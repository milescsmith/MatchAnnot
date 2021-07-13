# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2021-07-13

### Added
- `pyproject.toml`

### Changed
- reorganized code slightly to match module structure


## [20140728.01]
- Baselined. Output matches perl version.


## [20140804.01]
- Keep a single gene list per chr, rather than one per strand.
- Create AnnotationCursor class for walking annotation list.
- If not genes found, try the other strand.


## [20140805.01]
- Add multimap and opposite strand counts to summary lines.


## [20140821.01]
- Don't look for cluster ID if it isn't there.


## [20140903.01]
- Check that SAM file is sorted.
- Exons are numbered 1..N.


## [20140908.01] 
- Add # of exons to result line.
- Print cigar: line.


## [20140910.01] 
- Score=2 now means best 1-1 exon match.
- Print exons even if they don't match perfectly.


## [20140911.01] 
- Support for pickled annotation file.


## [20140922.01] 
- Revamp matchAnnot scoring system.


## [20140923.01] 
- Tidier version of cl: lines.


## [20140926.01] 
- Rename Clusters.py to ClusterReport.py.


## [20141015.01] 
- Accept alternate format GTF annotation file.
