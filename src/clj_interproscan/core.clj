(ns clj-interproscan.core
  (:require [clojure.data.xml :refer [parse]]
            [clojure.data.zip.xml :refer [xml1-> attr text xml-> tag=]]
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

(def tag-location?
  {:hmmer3-match :hmmer3-location
   :hmmer2-match :hmmer2-location
   :profilescan-match :profilescan-location
   :superfamilyhmmer3-match :superfamilyhmmer3-location
   :patternscan-match :patternscan-location
   :tmhmm-match :tmhmm-location
   :signalp-match :signalp-location
   :coils-match :coils-location
   :fingerprints-match :fingerprints-location
   :panther-match :panther-location
   :phobius-match :phobius-location
   :blastprodom-match :blastprodom-location
   :rpsblast-match :rpsblast-location})

(defmulti locations (fn [z lt]
                      (if (#{:patternscan-location :profilescan-location} lt)
                        :alignment
                        lt)))

(defmethod locations :default
  [z lt]
  (map :attrs (xml-> z :locations lt node)))

(defmethod locations :alignment
  [z lt]
  (map #(merge (:attrs (node %)) {:alignment (xml1-> % :alignment text)})
       (xml-> z :locations lt)))

(defn- signature-parser
  "Returns a map representing a hmmer match."
  [lt]
  #(let [sign (xml1-> % :signature)
         entry (xml1-> % :signature :entry)]
     (println lt)
     (merge (:attrs (node %))
            {:locations (locations % lt)}
            {:signature
             (merge (:attrs (node sign))
                    {:abstract (xml1-> sign :abstract text)}
                    {:comment (xml1-> sign :comment text)}
                    {:xrefs (map :attrs (xml-> sign :xref node))}
                    {:deprecated-acs (xml-> sign :deprecated-ac text)}
                    {:models (map :attrs (xml-> sign :models :model node))}
                    {:library (:attrs (xml1-> sign :signature-library-release node))}
                    {:entry
                     (if entry
                       (merge (:attrs (xml1-> entry node))
                              {:gos (map :attrs (xml-> entry :go-xref node))}
                              {:pathways (map :attrs (xml-> entry :pathway-xref node))}))})})))

(defn- parse-tag
  [z t]
  (map (signature-parser (tag-location? t)) (xml-> z :matches t)))

(defn hmmer-3-seq
  "Returns a lazy list of maps representing hmmer-3 matches."
  [zipper]
  (parse-tag zipper :hmmer3-match))

(defn hmmer-2-seq
  "Returns a lazy list of maps representing hmmer-2 matches."
  [zipper]
  (parse-tag zipper :hmmer2-match))

(defn profilescan-seq
  "Returns a lazy list of maps representing profilescan matches."
  [zipper]
  (parse-tag zipper :profilescan-match))

(defn superfamily-seq
  "Returns a lazy list of maps representing superfamily matches."
  [zipper]
  (parse-tag zipper :superfamilyhmmer3-match))

(defn patternscan-seq
  "Returns a lazy list of maps representing patternscan matches."
  [zipper]
  (parse-tag zipper :patternscan-match))

(defn tmhmm-seq
  "Returns a lazy list of maps representing tmhmm matches."
  [zipper]
  (parse-tag zipper :tmhmm-match))

(defn signalp-seq
  "Returns a lazy list of maps representing signalp matches."
  [zipper]
  (parse-tag zipper :signalp-match))

(defn coils-seq
  "Returns a lazy list of maps representing coils matches."
  [zipper]
  (parse-tag zipper :coils-match))

(defn fingerprints-seq
  "Returns a lazy list of maps representing fingerprint matches."
  [zipper]
  (parse-tag zipper :fingerprints-match))

(defn panther-seq
  "Returns a lazy list of maps representing panther matches."
  [zipper]
  (parse-tag zipper :panther-match))

(defn phobius-seq
  "Returns a lazy list of maps representing phobius matches."
  [zipper]
  (parse-tag zipper :phobius-match))

(defn blastprodom-seq
  "Returns a lazy list of maps representing blastprodom matches."
  [zipper]
  (parse-tag zipper :blastprodom-match))

(defn rpsblast-seq
  "Returns a lazy list of maps representing rpsblast matches."
  [zipper]
  (parse-tag zipper :rpsblast-match))

(defn any-seq
  "Returns a lazy list of maps representing all matches. The type of
  match is the value of :tag in the returned map."
  [zipper]
  (let [loc (first (xml-> zipper :matches node))]
    (->> (filter #((set (keys tag-location?)) (:tag %))
                 (:content loc))
         (map #(vector (xml-zip %) (:tag %)))
         (map (fn [[z t]]
                (assoc ((signature-parser (tag-location? t)) z) :tag t))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; go's
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn go-terms
  "Returns a lazy list of all GO entries from an interproscan zipper."
  [zipper]
  (->> (any-seq zipper)
       (mapcat #(get-in % [:signature :entry :gos]))))

(defn pathways
  "Returns a lazy list of all pathway entries from an interproscan zipper."
  [zipper]
  (->> (any-seq zipper)
       (mapcat #(get-in % [:signature :entry :pathways]))))

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
