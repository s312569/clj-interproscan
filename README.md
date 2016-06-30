# clj-interproscan

A Clojure library designed to for running and parsing Interproscan.

## Usage

To include add the following to your prjects.clj file:

```clojure
[clj-interproscan "0.1.0"]
```

Run interproscan on a file:

```clojure
user> (ips-file "/path/to/fasta/file" "/path/to/output/file")
Running: ["interproscan.sh" "-i" "/tmp/ips-input1467287656967-2885328136" "-o" "/path/to/output/file-1.xml" "-seqtype" "p" "-f" "XML" "-appl Pfam" "-dp" "-pa" "-iprlookup" "-goterms"]
user>
```

Sequences are run 10,000 at a time using pmap.

To parse open a buffered reader on the output file and call `ips-seq`
on the reader to prvide a lazy list of zippers representing
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

To access matches only `hmmer-3-seq` is definedat the moment which
returns a lazy sequence of maps containing match details. I'll add
more accessors as I need them but as 'ips-seq' returns a lazy list of
zippers it is easy to get at your information using
'clojure.data.zip.xml' and 'clojure.zip'.

```clojure
user> (with-open [r (io/reader tf)]
                         (->> (ips-seq r)
                              second
                              hmmer-3-seq
                              first))
{:evalue "8.5E-21", :score "73.6", :signature {:ac "PF07654", :desc "Immunoglobulin C1-set domain",
 :name "C1-set", :abstract nil, :comment nil, :xrefs (), :deprecated-acs (), :models ({:ac "PF07654",
 :desc "Immunoglobulin C1-set domain", :name "C1-set"}), :entry {:ac "IPR003597",
 :desc "Immunoglobulin C1-set", :name "Ig_C1-set", :type "DOMAIN", :gos (), :pathways ()}}}
user>
```
## License

Copyright © 2016 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
