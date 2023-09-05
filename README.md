# ROA

## Install
```sh
git clone git@github.com:dwpeng/roa.git
cd ./roa
make -j 4
```

## Index
```sh
./roa index ref.index ref1.fa ref2.fa
```

## Design
```sh
./roa design -i ref.index -q cDNA.fa
```

## Help message

```sh
# ./roa -h
ROA Template Designer.
Usage:
  ./roa <command> <options>
Commands:
  index         create index file
  design        design ROA template
```

```sh
# ./roa index -h
ROA Template Designer.
Usage:
  ./roa index output.index ref1.fa ref2.fa ...
Options:
  -h            show this help message
```

```sh
# ./roa design -h
ROA Template Designer.
Usage:
  ./roa design <options>
Example:
  ./roa design -i index.index -q query.fa -o template.fa -pairCheck 1
Options:
  -i <index>    index file path
  -q <query>    query file path
  -o <output>   output file path [template.fa]
  -homopolymer  homopolymer length [3]
  -minGC        min GC rate [0.45]
  -maxGC        max GC rate [0.55]
  -minTm        min melting temperature [52.4]
  -maxTm        max melting temperature [55.4]
  -avoidCGIn3   avoid CG in 3' end [1]
  -avoidTIn3    avoid T in 3' end [1]
  -pairCheck    check pair [0] maybe cost a long time
  -ncircle      number of circles [5]
  -h            show this help message
```
