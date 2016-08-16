# clj-interproscan

A Clojure library for running Interproscan and parsing the results.

## Usage

To include add the following to your projects.clj file:

```clojure
[clj-interproscan "0.2.0"]
```

To run interproscan on a fasta formatted file of protein sequences:

```clojure
user> (ips-file "/path/to/fasta/file" "/path/to/output/file")
Running: ["interproscan.sh" "-i" "/tmp/ips-input1467287656967-2885328136"
 "-o" "/path/to/output/file-1.xml" "-seqtype" "p" "-f" "XML" "-appl Pfam"
 "-dp" "-pa" "-iprlookup" "-goterms"]
user>
```

Sequences are run 10,000 at a time using pmap.

To parse open a buffered reader on the output file and call `ips-seq`
on the reader to provide a lazy list of zippers representing
interproscan protein entries:

```clojure
user> (with-open [r (io/reader tf)]
                         (->> (ips-seq r)
                              first
                              accession))
"P30447"
user> (with-open [r (io/reader tf)]
                         (->> (ips-seq r)
                              first
                              description))
"1A23_HUMAN HLA class I histocompatibility antigen, A-23 alpha chain OS=Homo sapiens GN=HLA-A PE=1 SV=1"
```

To access matches only `hmmer-3-seq` is defined at the moment which
returns a lazy sequence of maps containing match details. I'll add
more accessors as I need them but as 'ips-seq' returns a lazy list of
zippers it is easy to get at your information using
`clojure.data.zip.xml` and `clojure.zip`.

```clojure
user> (with-open [r (io/reader tf)]
                         (->> (ips-seq r)
                              first
                              hmmer-3-seq
                              first))
{:evalue "1.3E-13", :score "50.0", :signature {:ac "PF06623", :desc "MHC_I C-terminus",
 :name "MHC_I_C", :abstract nil, :comment nil, :xrefs (), :deprecated-acs (), :models
 ({:ac "PF06623", :desc "MHC_I C-terminus", :name "MHC_I_C"}), :entry {:ac "IPR010579",
 :desc "MHC class I, alpha chain, C-terminal", :name "MHC_I_a_C", :type "DOMAIN", :gos
 ({:category "BIOLOGICAL_PROCESS", :db "GO", :id "GO:0006955", :name "immune response"}
 {:category "CELLULAR_COMPONENT", :db "GO", :id "GO:0042612",
 :name "MHC class I protein complex"} {:category "CELLULAR_COMPONENT", :db "GO",
 :id "GO:0016020", :name "membrane"} {:category "BIOLOGICAL_PROCESS", :db "GO",
 :id "GO:0019882", :name "antigen processing and presentation"}), :pathways ()}}}
user>
```
## License

Copyright © 2016 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
