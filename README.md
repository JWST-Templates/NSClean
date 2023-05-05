# JWST NIRSpec Clean (NSClean)

27 March 2023 (Rev. 1.5)

## 1. Introduction

This software accompanies the draft PASP article that describes the algorithm. The current draft can be found in nsclean/doc/pasp_article.pdf.

NSClean is a post-processing tool for removing residual correlated read noise from NIRSpec images. These include "picture frame" noise and vertical banding. When not handled correctly, this correlated noise complicates calibration and can add spectral features that are not real. NSClean works by fitting and subtracting a background model that is constructed using areas of the NIRSpec focal plane that are either blanked-off, or thought to be usefully dark.

All background fitting is done in Fourier space. This facilitates: (1) using as many background pixels as possible and (2) fitting many degrees of freedom. The real power of NSClean  is reserved for observers who are able to specify background pixel masks. With few exceptions, masks are completely general. Observers who take the time to make thorough masks will be rewarded with very good correlated noise correction.

## 2. Implementation

NSClean is written in python-3. It was developed and tested using a scientific workstation having the following capabilities.

* 8x Intel(R) Xeon(R) CPUs E5-2637 v4 @ 3.50GHz CPUs
* 128 GB RAM (< 1GB is used by NSClean)
* 1x NVIDIA Quadro M4000 GPW w/ 8 GB RAM

In early versions of NSClean, using a GPU greatly accelerated computation. Subsequent algorithm improvements have made it so that NSClean runs about equally fast using CPUs and a GPU. The typical execution time on a mainstream workstation or laptop is a few seconds.

## 3. Theory

See the documentation folder

## 4. Examples

See the Notebooks folder.


## Modification history
March 2022, Bernard J. Rauscher, NASA Goddard Space Flight Center
* Conceived algorithm
* Revs. 1.0 - 1.5
* 1st release to JWST observers on 20 April 2023
