(ns clj-interproscan.core
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.java.io :as io]
            [clj-commons-exec :as exec]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clj-fasta.core :as fa]
            [me.raynes.fs :as fs]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; proteins
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn ips-seq
  "Returns a lazy list of 'protein' zippers from an Interproscan xml
  results file."
  [reader]
  (->> (:content (xml/parse reader))
       (filter (fn [x] (= :protein (:tag x))))
       (map zip/xml-zip)))

(defn accession
  [zipper]
  (zf/xml1-> zipper :xref (zf/attr :id)))

(defn description
  [zipper]
  (zf/xml1-> zipper :xref (zf/attr :desc)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; accessors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn hmmer-3-seq
  "Returns a map representing a hmmer-3 match."
  [zipper]
  (map #(let [sign (zf/xml1-> % :signature)
              entry (zf/xml1-> % :signature :entry)]
         (merge (:attrs (zip/node %))
                {:signature
                 (merge (:attrs (zip/node sign))
                        {:abstract (zf/xml1-> sign :abstract zf/text)}
                        {:comment (zf/xml1-> sign :comment zf/text)}
                        {:xrefs (map :attrs (zf/xml-> sign :xref zip/node))}
                        {:deprecated-acs (zf/xml-> sign :deprecated-ac zf/text)}
                        {:models (map :attrs (zf/xml-> sign :models :model zip/node))}
                        {:entry
                         (merge (:attrs (zf/xml1-> entry zip/node))
                                {:gos (map :attrs (zf/xml-> % :signature :entry :go-xref zip/node))}
                                {:pathways (map :attrs (zf/xml-> entry :pathway-xref zip/node))})})}))
       (zf/xml-> zipper :matches :hmmer3-match)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; running
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- ips-command
  [& {:keys [infile outfile appl precalc pathways lookup goterms seqtype]}]
  (vec
   (remove nil? (-> (list
                     "interproscan.sh" "-i" (str infile) "-o" (str outfile) "-seqtype" seqtype "-f" "XML"
                     (if appl (str "-appl " (apply str (interpose "," appl))))
                     (if (not precalc) "-dp")
                     (if pathways "-pa")
                     (if (or lookup goterms pathways) "-iprlookup")
                     (if goterms "-goterms"))))))

(defn- run-ips
  [& {:keys [infile outfile appl precalc pathways lookup goterms seqtype]}]
  (try
    (let [command (ips-command :infile infile
                               :outfile outfile
                               :appl appl
                               :precalc precalc
                               :seqtype seqtype
                               :lookup lookup
                               :goterms goterms
                               :pathways pathways)
          ips (do
                (println (str "Running: " command))
                @(exec/sh command))]
      (if (= 0 (:exit ips))
        outfile
        (do (println ips)
            (throw (Throwable. (str "Interproscan error: " (:err ips)))))))
    (catch Exception e
      (fs/delete outfile)
      (throw e))))

(defn ips
  "Runs interproscan on a collection of fasta formatted protein
  sequences (see clj-fasta). Specify analyses using the :appl keyword,
  default is \"Pfam\" only. To run all analyses set :appl to
  nil. Splits sequences into lots of 10,000 and runs interproscan
  concurrently on each group using pmap."
  [coll outfile {:keys [appl lookup goterms precalc pathways seqtype]
                 :or {appl '("Pfam") lookup true goterms true precalc false pathways true seqtype "p"}}]
  (let [c (atom 0)]
    (pmap 
     #(let [i (fa/fasta->file % (fs/temp-file "ips-input") :append false)]
        (try
          (run-ips :infile (fs/absolute i)
                   :outfile (str (fs/absolute outfile) "-" (swap! c inc) ".xml")
                   :appl appl
                   :precalc precalc
                   :pathways pathways
                   :seqtype seqtype
                   :lookup lookup
                   :goterms goterms)
          (finally (fs/delete i))))
     (partition-all 10000 coll))))

(defn ips-file
  "Runs interproscan on a file of fasta formatted protein
  sequences. Specify analyses using the :appl keyword, default is
  \"Pfam\" only. To run all analyses set :appl to nil. Splits
  sequences into lots of 10,000 and runs interproscan concurrently on
  each group using pmap."
  ([file outfile] (ips-file file outfile {}))
  ([file outfile {:keys [appl lookup goterms precalc pathways seqtype]
                  :or {appl '("Pfam") lookup true goterms true precalc false pathways true seqtype "p"}
                  :as m}]
   (with-open [r (io/reader file)]
     (ips (fa/fasta-seq r) outfile m))))
