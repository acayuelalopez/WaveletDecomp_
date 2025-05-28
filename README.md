# WaveletDecomp Plugin for ImageJ

**WaveletDecomp** is a plugin for ImageJ designed to facilitate wavelet decomposition and reconstruction of images. It provides a graphical interface for selecting, processing, and analyzing image data using wavelet transforms.

## Main Features

- Load and preview multiple images from a directory.
- Perform forward and inverse wavelet transforms.
- Automatic and manual modes for wavelet analysis.
- Calculation of wavelet coefficients and reconstruction of images.
- Export results in CSV format.
- Save processed images as TIFF files.

## Requirements

- Java 8 or higher.
- ImageJ (preferably the Fiji distribution).
- Apache Commons Lang library for array manipulation.

## Installation

1. Compile or download the JAR file of the plugin.
2. Place it in the `plugins` folder of ImageJ/Fiji.
3. Restart ImageJ.
4. Access the plugin via `Plugins > WaveletDecomp`.

## Workflow

1. Select a directory of images.
2. Define a region of interest (ROI) with the rectangle tool.
3. Perform forward wavelet transform to decompose the image.
4. Perform inverse wavelet transform to reconstruct the image.
5. Export the results and save the processed images.

## Outputs

- CSV file with quantitative data per image.
- TIFF images of the masks and generated montages.
- Individual and total results tables.

## Application

This plugin is useful for researchers in image processing and analysis, enabling semi-automated wavelet decomposition and reconstruction of images.

## License

Distributed under the MIT License.

