# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


Version cheatsheet: <MAJOR>.<MINOR>.<PATCH>

Patches should not contain new functionality, only bugfixes.
Minor version should not contain breaking changes.

Types of changes:

* Added
* Changed
* Deprecated
* Removed
* Fixed
* Security

## [Unreleased]

## [0.2.0] - 2024-05-23

### Added

- Adding a option to simulate a simple recombination with a specific number of cross-overs per chromosome (some organisms have tight regulations, e.g. C. elegans)
- A test to control that cross-overs are returned sorted per segment

### Fixed

- Fixed order of crossovers

### Removed

- The option to not export in dwgsim format as it was not implemented, and likely wont

## [0.1.0] - 2024-02-18

First experimental version, unfinished.
