# Changelog

## [0.3.1](https://github.com/hydra-genetics/reports/compare/v0.3.0...v0.3.1) (2024-01-12)


### Bug Fixes

* address transition bugs when toggling callers ([#57](https://github.com/hydra-genetics/reports/issues/57)) ([1154903](https://github.com/hydra-genetics/reports/commit/1154903f7cc86b8f575eed48f3683054fa33a644))

## [0.3.0](https://github.com/hydra-genetics/reports/compare/v0.2.0...v0.3.0) (2024-01-08)

This release contains quite a few updates. Since the module currently only contains the CNV report, all changes are related to this. Perhaps the most significant of these changes is that the report template files are now accessed from the reports module itself, so it is no longer necessary to include the template in the parent workflow. Information on how to customize the template is included in [the documentation](https://hydra-genetics-reports.readthedocs.io/en/latest/reports/).

Previously it was not possible to view the resulting reports without a working internet connection, but now all external dependencies have been integrated into the report. The only requirement is a reasonably modern web browser.

An attempt to improve the performance of the plotting in the CNV report has been done. Now, if there are more than a certain number of points being displayed per chromosome (currently 300), the data is binned to display the median &pm; 1 standard deviation. For the VAF-plots it makes an attempt to detect whether the distribution within a bin is bimodal, in which case it plots two values: one representing values &gt; 0.5 and one representing values &leq; 0.5. Zooming in so that the number of points are fewer than this limit will result in all points being shown, both in the case of the log<sub>2</sub> ratio and the VAF.

### Features

* add template version check ([#46](https://github.com/hydra-genetics/reports/issues/46)) ([dd780f0](https://github.com/hydra-genetics/reports/commit/dd780f0efc0236584b8fdee91bbd0853ab0ee677))
* embed external dependencies ([#51](https://github.com/hydra-genetics/reports/issues/51)) ([55a05bd](https://github.com/hydra-genetics/reports/commit/55a05bd77a690beb065b7204db881f6baaf5adaa))
* improve plotting performance ([#42](https://github.com/hydra-genetics/reports/issues/42)) ([87f9542](https://github.com/hydra-genetics/reports/commit/87f9542539c7589de9336a7a275aa025e92fff4f))
* use template directly from the reports repo ([#48](https://github.com/hydra-genetics/reports/issues/48)) ([7848269](https://github.com/hydra-genetics/reports/commit/7848269f43444fc55091139b37ebf34568417e32))


### Bug Fixes

* **cnv:** bug in checking input keys ([#49](https://github.com/hydra-genetics/reports/issues/49)) ([af34098](https://github.com/hydra-genetics/reports/commit/af340981e520e6e95f9649622d6b8e1e94ebce9a))
* handle huge csv files ([#43](https://github.com/hydra-genetics/reports/issues/43)) ([73c9043](https://github.com/hydra-genetics/reports/commit/73c904387f0c40c54c62dd2676b58534015644ab))
* optimise transitions in plots ([#45](https://github.com/hydra-genetics/reports/issues/45)) ([651c8c0](https://github.com/hydra-genetics/reports/commit/651c8c0d5fee3506707237374855e0123fceca89))

## [0.2.0](https://github.com/hydra-genetics/reports/compare/v0.1.0...v0.2.0) (2023-09-21)


### Features

* add cytoband information to chromosome plot ([#35](https://github.com/hydra-genetics/reports/issues/35)) ([620cae5](https://github.com/hydra-genetics/reports/commit/620cae58ac35abcb102487039cebe3bb2ff170ee))
* limit maximum zoom in chromosome plot ([#31](https://github.com/hydra-genetics/reports/issues/31)) ([9aed5a1](https://github.com/hydra-genetics/reports/commit/9aed5a1a4a99e62428872a710828c5d3e87f2aab))


### Bug Fixes

* add tc-file to cnv_html_report input to avoid that it is removed ([#39](https://github.com/hydra-genetics/reports/issues/39)) ([b891111](https://github.com/hydra-genetics/reports/commit/b891111a835ffdd201ce405151014f11ff545ede))
* issue with zooming into region with no data ([#30](https://github.com/hydra-genetics/reports/issues/30)) ([1107449](https://github.com/hydra-genetics/reports/commit/11074497e57652c3d2116260ae0642f3989fb88b))
* update chromosome lengths properly ([#32](https://github.com/hydra-genetics/reports/issues/32)) ([d8b41c8](https://github.com/hydra-genetics/reports/commit/d8b41c80dddf9817f73849908c66bf005d0cc600))
* wrap gene names in table ([#37](https://github.com/hydra-genetics/reports/issues/37)) ([11d6c3a](https://github.com/hydra-genetics/reports/commit/11d6c3a9dd1e22352f0c7021b352e8049ed2c175))

## 0.1.0 (2023-06-02)


### Features

* add CNV VCF support ([2312e8f](https://github.com/hydra-genetics/reports/commit/2312e8f3d1c4ba287fced98467cccebcd42fe365))
* add support for a germline vcf ([abd7ad3](https://github.com/hydra-genetics/reports/commit/abd7ad34d64a45196605bb62080ffe1c572c017a))
* supply bed file with annotations for chromosome plot ([91daf6b](https://github.com/hydra-genetics/reports/commit/91daf6bc9b199b5e38c8d7c80d2cd27215f1d9c0))


### Bug Fixes

* correct errors introduced by prettier ([a4a3c7e](https://github.com/hydra-genetics/reports/commit/a4a3c7e3402d2c0f9b1c9b3ecd158592427c4c8c))
* crash when CNVs had more than one gene ([e3166fc](https://github.com/hydra-genetics/reports/commit/e3166fcf2246d6fe813ea484dde4b8c42b2e74ab))
* don't generate copy rules if used as module ([88e7941](https://github.com/hydra-genetics/reports/commit/88e79419600e66c4c1ba20f003b61fd7e5ede6d0))
* faulty coordinates for annotations with same name ([bde1576](https://github.com/hydra-genetics/reports/commit/bde15767469b367857b3cbd5a171119e6a1c3363))
* incorrect ratio file for GATK ([41b44eb](https://github.com/hydra-genetics/reports/commit/41b44eba3548e3059f5768e78b65af35d3d87b22))
* issue with updates when changing chromosomes ([ea04e73](https://github.com/hydra-genetics/reports/commit/ea04e738e04469f1c4edbed02c4a81a2a673ab89))


### Miscellaneous Chores

* release 0.1.0 ([#16](https://github.com/hydra-genetics/reports/issues/16)) ([e02f7e2](https://github.com/hydra-genetics/reports/commit/e02f7e29e77b2b710fee3a6926ce57b05fd7590f))
