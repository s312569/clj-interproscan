(ns clj-interproscan.core
  (:require [clojure.data.xml :refer [parse]]
            [clojure.data.zip.xml :refer [xml1-> attr text xml->]]
            [clojure.java.io :refer [reader]]
            [clj-commons-exec :refer [sh]]
            [clojure.zip :refer [xml-zip node]]
            [clojure.string :refer []]
            [clj-fasta.core :refer [fasta->file fasta-seq]]
            [biodb.core :refer [table-spec prep-sequences freeze thaw restore-sequence]]
            [me.raynes.fs :refer [delete temp-file absolute]]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; proteins
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn ips-seq
  "Returns a lazy list of 'protein' zippers from an Interproscan xml
  results file."
  [reader]
  (->> (:content (parse reader))
       (filter (fn [x] (= :protein (:tag x))))
       (map xml-zip)))

(defn accession
  [zipper]
  (xml1-> zipper :xref (attr :id)))

(defn description
  [zipper]
  (xml1-> zipper :xref (attr :desc)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; accessors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn hmmer-3-seq
  "Returns a map representing a hmmer-3 match."
  [zipper]
  (map #(let [sign (xml1-> % :signature)
             entry (xml1-> % :signature :entry)]
         (merge (:attrs (node %))
                {:signature
                 (merge (:attrs (node sign))
                        {:abstract (xml1-> sign :abstract text)}
                        {:comment (xml1-> sign :comment text)}
                        {:xrefs (map :attrs (xml-> sign :xref node))}
                        {:deprecated-acs (xml-> sign :deprecated-ac text)}
                        {:models (map :attrs (xml-> sign :models :model node))}
                        {:entry
                         (if entry
                           (merge (:attrs (xml1-> entry node))
                                  {:gos (map :attrs (xml-> entry :go-xref node))}
                                  {:pathways (map :attrs (xml-> entry :pathway-xref node))}))})}))
       (xml-> zipper :matches :hmmer3-match)))

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
          ips (do (println (str "Running: " command)) @(sh command))]
      (if (= 0 (:exit ips))
        outfile
        (throw (Exception. (str "Interproscan error: " (:err ips))))))))

(defn ips
  "Runs interproscan on a collection of fasta formatted protein
  sequences (see clj-fasta). Specify analyses using the :appl keyword,
  default is \"Pfam\" only. To run all analyses set :appl to
  nil. Splits sequences into lots of 10,000 and runs interproscan
  on each group."
  [coll outfile {:keys [appl lookup goterms precalc pathways seqtype]
                 :or {appl '("Pfam") lookup true goterms true precalc false pathways true seqtype "p"}}]
  (let [c (atom 0)
        fl (atom [])]
    (try
      (doall
       (map 
        #(let [i (fasta->file % (temp-file "ips-input") :append false)
               o (str (absolute outfile) "-" (swap! c inc) ".xml")]
           (swap! fl conj o)
           (try
             (run-ips :infile (absolute i)
                      :outfile o
                      :appl appl
                      :precalc precalc
                      :pathways pathways
                      :seqtype seqtype
                      :lookup lookup
                      :goterms goterms)
             (finally (delete i))))
        (partition-all 10000 coll)))
      (catch Exception e
        (doseq [f @fl] (delete f))
        (throw e)))))

(defn ips-file
  "Runs interproscan on a file of fasta formatted protein
  sequences. Specify analyses using the :appl keyword, default is
  \"Pfam\" only. To run all analyses set :appl to nil. Splits
  sequences into lots of 10,000 and runs interproscan on each group."
  ([file outfile] (ips-file file outfile {}))
  ([file outfile {:keys [appl lookup goterms precalc pathways seqtype]
                  :or {appl '("Pfam") lookup true goterms true precalc false pathways true seqtype "p"}
                  :as m}]
   (with-open [r (reader file)]
     (ips (fasta-seq r) outfile m))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; integration with biodb
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod table-spec :ips
  [q]
  (vector [:accession :text "PRIMARY KEY"]
          [:src :binary "NOT NULL"]))

(defmethod prep-sequences :ips
  [q]
  (->> (:coll q)
       (map #(hash-map :accession (accession %) :src (freeze (node %))))))

(defmethod restore-sequence :ips
  [q]
  (xml-zip (thaw (:src (dissoc q :type)))))
