import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.ImageStack;
import ij.CompositeImage;
import ij.measure.Calibration;
import ij.plugin.filter.FFTFilter;
import ij.IJ;

import java.awt.Point;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.Raster;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;

import ij.ImagePlus;

class ImageData {
	private final boolean DEBUG = false;
	protected ImagePlus imageOrig;
	protected ImagePlus imageWave;
	protected ImagePlus imageModif;
	protected ImagePlus imageWaveNoStretch;
	private int type;
	private int width;
	private int height;
	private int nChannels;
	private int nSlices;
	private int imageSize;
	private double[][][] imageData;
	private double[][][] transformedData;
	private int scale;
	private double lowCoeffA;
	private double lowCoeffD;
	private double highCoeffA;
	private double highCoeffD;
	private double[] minInputVals;
	private double[] maxInputVals;
	private double minOutputVal;
	private double maxOutputVal;
	private WaveletFilters waveletFilter;
	private BufferedWriter bw;
	double[] detailDef;
	double[] approx, detail, approxInitial, detailVert, detailDiag, detailHoriz;
	Calibration cal;
	byte[] bytes;

	protected int getNSlices() {
		return this.nSlices;
	}

	protected int getScale() {
		return this.scale;
	}

	protected void setScale(final int scale) {
		this.scale = scale;
	}

	protected int getType() {
		return this.type;
	}

	protected int getHeight() {
		return this.height;
	}

	protected int getWidth() {
		return this.width;
	}

	protected ImageData(final ImagePlus imp) {
		this.imageOrig = imp;
		this.type = this.imageOrig.getType();
		this.width = this.imageOrig.getWidth();
		this.height = this.imageOrig.getHeight();
		this.nSlices = this.imageOrig.getStackSize();
		this.nChannels = this.imageOrig.getNChannels();
		this.imageSize = this.width * this.height;
		this.imageData = new double[this.nSlices][this.height][this.width];
		this.transformedData = new double[this.nSlices][this.height][this.width];
		this.minInputVals = new double[this.nSlices];
		this.maxInputVals = new double[this.nSlices];
		this.waveletFilter = new WaveletFilters();
		this.imageWave = this.imageOrig.duplicate();
		this.imageModif = this.imageOrig.duplicate();
		cal = this.imageOrig.getCalibration();
		(this.imageWaveNoStretch = IJ.createImage("WT-NoStretch", this.width, this.height, this.nSlices, 32))
				.setCalibration(cal);

//		this.bitmapToData();
//		this.dataToBitmap();

	}

	protected void bitmapToData(final int z) {

		this.image2DToData(z, (byte[]) this.imageOrig.getStack().getProcessor(z).getPixels());

	}

	private void prepareStretchOI(final int z) {
		this.minOutputVal = Double.MAX_VALUE;
		this.maxOutputVal = Double.MIN_VALUE;
		for (int y = 0; y < this.height; ++y) {
			for (int x = 0; x < this.width; ++x) {
				final double outputVal = this.imageData[z - 1][y][x];
				if (this.minOutputVal > outputVal) {
					this.minOutputVal = outputVal;
				}
				if (this.maxOutputVal < outputVal) {
					this.maxOutputVal = outputVal;
				}
			}
		}
	}

	private double stretchOI(final int z, final double x) {
		final double outMax = this.maxInputVals[z - 1];
		final double outMin = this.minInputVals[z - 1];
		final double inMax = this.maxOutputVal;
		final double inMin = this.minOutputVal;
		double maxGreyVal = 255.0;
		if (this.type == 2) {
			final double scale = (outMax - outMin) / (inMax - inMin);
			return scale * (x - inMin) + outMin;
		}
		if (this.type == 0) {
			maxGreyVal = 255.0;
		}
		if (this.type == 1) {
			maxGreyVal = 65535.0;
		}
		if (inMin != inMax) {
			final double scale = (outMax - outMin) / (inMax - inMin);
			return scale * (x - inMin) + outMin;
		}
		if (x > maxGreyVal) {
			return maxGreyVal;
		}
		if (x >= 0.0 && x <= maxGreyVal) {
			return x;
		}
		return 0.0;
	}

	private void image2DToData(final int z, final byte[] pixels) {
		double minVal = Double.MAX_VALUE;
		double maxVal = Double.MIN_VALUE;
		for (int y = 0; y < this.height; ++y) {
			for (int x = 0; x < this.width; ++x) {
				final int index = x + y * this.width;
				final byte pixelByte = pixels[index];
				final double pixelDouble = pixelByte & 0xFF;
				this.imageData[z - 1][y][x] = pixelDouble;
				if (minVal > pixelDouble)
					minVal = pixelDouble;

				if (maxVal < pixelDouble)
					maxVal = pixelDouble;

			}
		}
		this.minInputVals[z - 1] = minVal;
		this.maxInputVals[z - 1] = maxVal;
	}

	protected void dataToBitmap(final int z) {

		this.dataToImage2D(z, (byte[]) this.imageModif.getStack().getProcessor(z).getPixels());

	}

	private void dataToImage2D(final int z, final byte[] pixels) {

		this.prepareStretchOI(z);
		for (int y = 0; y < this.height; ++y) {
			for (int x = 0; x < this.width; ++x) {
				final int index = x + y * this.width;
				pixels[index] = (byte) this.stretchOI(z, this.imageData[z - 1][y][x]);

			}
			// new ImagePlus(z+"", new ByteProcessor(this.width, this.height,
			// pixels)).show();
		}
	}

	private void dataToWave(final int z) {
		int row = 0;
		int col = 0;
		final int borderX = this.width / (int) Math.pow(2.0, this.scale);
		final int borderY = this.height / (int) Math.pow(2.0, this.scale);
		// if (Wavelet_Denoise.isShownWT()) {
		final double[][] tempDataStretched = new double[this.height][this.width];
		final double maxVal = 255.0;
		for (int i = 0; i < this.imageSize; ++i) {
			final double intensityStretched = Math.min(maxVal, Math.max(0.0,
					this.stretchWC(this.transformedData[z - 1][row][col], row < borderY && col < borderX)));
			tempDataStretched[row][col] = intensityStretched;
			if (++col == this.width) {
				col = 0;
				++row;
			}
		}

		this.dataToImage2D(tempDataStretched, (byte[]) this.imageWave.getStack().getProcessor(z).getPixels());

		// }
		// if (Wavelet_Denoise.isShownNoStretchWT()) {
		final double[][] tempDataNoStretch = new double[this.height][this.width];
		row = 0;
		col = 0;
		for (int j = 0; j < this.imageSize; ++j) {
			tempDataNoStretch[row][col] = this.transformedData[z - 1][row][col];
			if (++col == this.width) {
				col = 0;
				++row;
			}
		}
		this.dataToImage2D(tempDataNoStretch, (float[]) this.imageWaveNoStretch.getStack().getProcessor(z).getPixels());

	}

	private void dataToImage2D(final double[][] slice, byte[] pixels) {
		for (int y = 0; y < this.height; ++y) {
			for (int x = 0; x < this.width; ++x) {
				final int index = x + y * this.width;
				pixels[index] = (byte) slice[y][x];
				// new ImagePlus("", new ByteProcessor(this.width, this.height, pixels)).show();
			}
		}
	}

	private void dataToImage2D(final double[][] slice, final float[] pixels) {
		for (int y = 0; y < this.height; ++y) {
			for (int x = 0; x < this.width; ++x) {
				final int index = x + y * this.width;
				pixels[index] = (float) slice[y][x];
			}
		}
	}

	protected int changeScale(final int scale) {
		final int maxScale = (int) (Math.log10((this.width < this.height) ? this.width : this.height)
				/ Math.log10(2.0));
		if (scale < 1 || scale > maxScale) {
			return 0;
		}
		this.setScale(scale);
		return maxScale;
	}

	protected void setWaveletFilter(final String filterMenuName) {
		Label_2610: {
			switch (filterMenuName) {
			case "Haar 1": {
				this.waveletFilter.setFilter("haar1.flt");
				break Label_2610;
			}
			}
			this.waveletFilter.setFilter("haar1.flt");
		}
		this.transposeFilters();
	}

	private void transposeFilters() {
		final int len = this.waveletFilter.dlf.length;
		final double[] dlfT = new double[len];
		final double[] dhfT = new double[len];
		final double[] rlfT = new double[len];
		final double[] rhfT = new double[len];
		for (int i = 0; i < this.waveletFilter.dlf.length; ++i) {
			dlfT[i] = this.waveletFilter.dlf[len - i - 1];
			dhfT[i] = this.waveletFilter.dhf[len - i - 1];
			rlfT[i] = this.waveletFilter.rlf[len - i - 1];
			rhfT[i] = this.waveletFilter.rhf[len - i - 1];
		}
		this.waveletFilter.dlf = dlfT;
		this.waveletFilter.dhf = dhfT;
		this.waveletFilter.rlf = rlfT;
		this.waveletFilter.rhf = rhfT;
	}

	private void prepareStretchWC(final int z) {
		// if (Wavelet_Denoise.isShownWT()) {
		final int borderX = this.width / (int) Math.pow(2.0, this.scale);
		final int borderY = this.height / (int) Math.pow(2.0, this.scale);
		final double[] approx = new double[borderX * borderY];
		final double[] detail = new double[this.imageSize - approx.length];
		detailDef = new double[this.imageSize - approx.length];
		int indexA = 0;
		int indexD = 0;
		for (int y = 0; y < this.height; ++y) {
			for (int x = 0; x < this.width; ++x) {
				if (y < borderY && x < borderX) {
					approx[indexA] = this.transformedData[z - 1][y][x];
					++indexA;
				} else {
					detailDef[indexD] = this.transformedData[z - 1][y][x];
					detail[indexD] = this.transformedData[z - 1][y][x];
					++indexD;
				}
			}
		}
		final int percentile = 100;
		Arrays.sort(approx);
		this.lowCoeffA = approx[approx.length / percentile];
		this.highCoeffA = approx[approx.length - 1 - approx.length / percentile];
		Arrays.sort(detail);
		this.lowCoeffD = detail[detail.length / percentile];
		this.highCoeffD = detail[detail.length - 1 - detail.length / percentile];
		// }
	}

	private double stretchWC(final double x, final boolean isApprox) {
		final double outMax = 255.0;
		final double outMin = 0.0;
		final double inMax = isApprox ? this.highCoeffA : this.highCoeffD;
		final double inMin = isApprox ? this.lowCoeffA : this.lowCoeffD;
		if (inMin != inMax) {
			final double scale = (outMax - outMin) / (inMax - inMin);
			return scale * (Math.min(inMax, Math.max(inMin, x)) - inMin) + outMin;
		}
		if (x > 0.0) {
			return 255.0;
		}
		return 0.0;
	}

	private void FWT(final double[] data) {
		final int dataLen = data.length;
		final int kernelLen = this.waveletFilter.dlf.length;
		final int mid = dataLen / 2;
		final int midKernel = kernelLen / 2;
		final int pad = dataLen * 30;
		approx = new double[dataLen];
		detail = new double[dataLen];

		for (int i = 0; i < dataLen; ++i) {
			double sumL = 0.0;
			double sumH = 0.0;
			for (int j = 0; j < kernelLen; ++j) {

				sumL += data[(i - midKernel + j + pad) % dataLen] * this.waveletFilter.dlf[j];
				sumH += data[(i - midKernel + j + pad) % dataLen] * this.waveletFilter.dhf[j];
			}

			approx[i] = sumL;
			detail[i] = sumH;

		}

		for (int k = 0; k < mid; ++k) {
			final int l = k * 2;
			data[k] = approx[l];
			data[mid + k] = detail[l];
		}

	}

//	private void FWTNoFilter(final double[] data) {
//		final int dataLen = data.length;
//		final int kernelLen = this.waveletFilter.dlf.length;
//		final int mid = dataLen / 2;
//		final int midKernel = kernelLen / 2;
//		final int pad = dataLen * 30;
//		approx = new double[dataLen];
//		detail = new double[dataLen];
//
//		for (int i = 0; i < dataLen; ++i) {
//			double sumL = 0.0;
//			double sumH = 0.0;
//			for (int j = 0; j < kernelLen; ++j) {
//
//				sumL += data[(i - midKernel + j + pad) % dataLen] * this.waveletFilter.dlf[j];
//				sumH += data[(i - midKernel + j + pad) % dataLen] * this.waveletFilter.dhf[j];
//			}
//
//			approx[i] = sumL;
//			detail[i] = sumH;
//
//		}
//
//		for (int k = 0; k < mid; ++k) {
//			final int l = k * 2;
//			data[k] = approx[l];
//			data[mid + k] = detail[l];
//		}
//
//	}

	private void FWT2D(final int z) {
		final double[][] tempData = new double[this.height][this.width];

		for (int y = 0; y < this.height; ++y) {
			for (int x = 0; x < this.width; ++x) {
				tempData[y][x] = this.imageData[z - 1][y][x];
			}
		}

		for (int k = 0; k < this.scale; ++k) {
			final int lev = 1 << k;
			final int levCols = this.height / lev;
			final int levRows = this.width / lev;

			final double[] row = new double[levCols];

			for (int x2 = 0; x2 < levRows; ++x2) {
				for (int y2 = 0; y2 < row.length; ++y2) {
					row[y2] = tempData[y2][x2];
					// IJ.log("row: "+tempData[y2][x2]);
				}

				this.FWT(row);
				for (int y2 = 0; y2 < row.length; ++y2) {
					tempData[y2][x2] = row[y2];
					// IJ.log("row2: "+tempData[y2][x2]);

				}
			}
			final double[] col = new double[levRows];
			// byte[] pixels = new byte[levRows * levCols];
			byte[] pixels0 = new byte[this.imageSize];

			for (int y3 = 0; y3 < levCols; ++y3) {

				for (int x3 = 0; x3 < col.length; ++x3) {
					col[x3] = tempData[y3][x3];
					// IJ.log("col: "+tempData[y3][x3]);
				}
				this.FWT(col);

				for (int x3 = 0; x3 < col.length; ++x3) {

					// final int index = x3 + y3 * levRows;
					final int index0 = x3 + y3 * this.width;
					//tempData[y3][x3] = col[x3];
					// pixels[index] = (byte) col[x3];

					pixels0[index0] = Double.valueOf(col[x3]).byteValue();

				}

				// IJ.log(this.minInputVals[z - 1] + "----" + this.maxInputVals[z - 1]);
				ImagePlus impToFilter = new ImagePlus("", new ByteProcessor(this.width, this.height, pixels0));
				impToFilter.setRoi(levRows / 2, 0, levRows / 2, levCols / 2);
				// impToFilter.show();
				IJ.run(impToFilter, "Bandpass Filter...",
						String.format("filter_large=%d filter_small=1 suppress=Vertical tolerance=95",
								(int) ((impToFilter.getWidth() * 10) / 100)));
				impToFilter.killRoi();
				bytes = (byte[]) impToFilter.getProcessor().getPixels();

				for (int x3 = 0; x3 < col.length; ++x3) {

					tempData[y3][x3] = (Double.valueOf(bytes[x3]).doubleValue() / 255)
							* (this.maxInputVals[z - 1] - this.minInputVals[z - 1]) + this.minInputVals[z - 1];

				}

			}

			ImagePlus impApprox, impDetailVert, impDetailHoriz, impDetailDiag;
			ImagePlus impTotal = new ImagePlus("Total-", new ByteProcessor(this.width, this.height, pixels0));
			impTotal.setRoi(0, 0, levRows, levCols);
			ImagePlus impTotal2 = impTotal.crop();

			ImagePlus impToRoi = impTotal2.duplicate();

			ImagePlus impApproxCrop = impToRoi.duplicate();
			impApproxCrop.setRoi(0, 0, levRows / 2, levCols / 2);
			impApprox = impApproxCrop.crop();

			ImagePlus impDetailVertCrop = impToRoi.duplicate();
			impDetailVertCrop.setRoi(0, levCols / 2, levRows / 2, levCols / 2);
			impDetailVert = impDetailVertCrop.crop();

			ImagePlus impDetailHorizCrop = impToRoi.duplicate();
			impDetailHorizCrop.setRoi(levRows / 2, 0, levRows / 2, levCols / 2);
			impDetailHoriz = impDetailHorizCrop.crop();
			// IJ.run(impDetailHoriz, "Bandpass Filter...",
			// String.format("filter_large=%d filter_small=1 suppress=Vertical
			// tolerance=95",
			// (int) ((impDetailHoriz.getWidth() * 10) / 100)));
			ImagePlus impDetailDiagCrop = impToRoi.duplicate();
			impDetailDiagCrop.setRoi(levRows / 2, levCols / 2, levRows / 2, levCols / 2);
			impDetailDiag = impDetailDiagCrop.crop();

			IJ.saveAsTiff(impTotal2, "/home/anaacayuela/Ana_pruebas_imageJ/javi_conesa/jonathan/images/result_prueba"
					+ File.separator + "Total-" + k + "-" + z);

			IJ.saveAsTiff(impApprox, "/home/anaacayuela/Ana_pruebas_imageJ/javi_conesa/jonathan/images/result_prueba"
					+ File.separator + "Approx-" + k + "-" + z);
			IJ.saveAsTiff(impDetailVert,
					"/home/anaacayuela/Ana_pruebas_imageJ/javi_conesa/jonathan/images/result_prueba" + File.separator
							+ "Vert-" + k + "-" + z);
			IJ.saveAsTiff(impDetailHoriz,
					"/home/anaacayuela/Ana_pruebas_imageJ/javi_conesa/jonathan/images/result_prueba" + File.separator
							+ "Horiz-" + k + "-" + z);
			IJ.saveAsTiff(impDetailDiag,
					"/home/anaacayuela/Ana_pruebas_imageJ/javi_conesa/jonathan/images/result_prueba" + File.separator
							+ "Diag-" + k + "-" + z);

		}
		for (int y4 = 0; y4 < this.height; ++y4) {
			for (int x4 = 0; x4 < this.width; ++x4) {

				this.transformedData[z - 1][y4][x4] = tempData[y4][x4];

			}
		}
	}

//	private void FWT2DToFilter(final int z, ImagePlus impToFilter) {
//		final double[][] tempData = new double[this.height][this.width];
//		bytes = (byte[]) impToFilter.getProcessor().getPixels();
//		double minVal = Double.MAX_VALUE;
//		double maxVal = Double.MIN_VALUE;
//		for (int y = 0; y < this.height; ++y) {
//			for (int x = 0; x < this.width; ++x) {
//				final int index = x + y * this.width;
//				final byte pixelByte = bytes[index];
//				final double pixelDouble = pixelByte & 0xFF;
//				tempData[y][x] = pixelDouble;
//				// + "-------valuesProcesso: "
////							+ +Double.valueOf(tempData[y3][x3]).byteValue() + "------"
////							+ Byte.valueOf(Double.valueOf(tempData[y3][x3]).byteValue()).doubleValue()
////							+ "----- se guarda??:   " + Double.valueOf(((int) tempData[y3][x3])).doubleValue() + "----"
////							+ (((int) Double.valueOf(tempData[y3][x3]).doubleValue()) - 256));
//				if (minVal > pixelDouble)
//					minVal = pixelDouble;
//
//				if (maxVal < pixelDouble)
//					maxVal = pixelDouble;
//
//			}
//		}
//		this.minInputVals[z - 1] = minVal;
//		this.maxInputVals[z - 1] = maxVal;
//		byte[] bytes = new byte[this.imageSize];
//		for (int k = 0; k < this.scale; ++k) {
//			final int lev = 1 << k;
//			final int levCols = this.height / lev;
//			final int levRows = this.width / lev;
//			final double[] row = new double[levCols];
//	
////			for (int x2 = 0; x2 < levRows; ++x2) {
////				for (int y2 = 0; y2 < row.length; ++y2) {
////					row[y2] = tempData[y2][x2];
////				}
////				// this.FWTNoFilter(row);
////				for (int y2 = 0; y2 < row.length; ++y2) {
////					tempData[y2][x2] = row[y2];
////					final int index0 = y2 + x2 * this.width;
////					bytes[index0] = (byte) row[y2];
////				}
////			}
////			final double[] col = new double[levRows];
////		
////			for (int y3 = 0; y3 < levCols; ++y3) {
////				for (int x3 = 0; x3 < col.length; ++x3) {
////					col[x3] = tempData[y3][x3];
////				}
////				//this.FWTNoFilter(col);
////				for (int x3 = 0; x3 < col.length; ++x3) {
////					tempData[y3][x3] = col[x3];
////					final int index0 = x3 + y3 * this.width;
////					bytes[index0] = (byte)col[x3];
////					//IJ.log("tempData: " + Double.valueOf(tempData[y3][x3]).doubleValue());
////				}
////			}
//			for (int y4 = 0; y4 < this.height; ++y4) {
//				for (int x4 = 0; x4 < this.width; ++x4) {
//					final int index0 = x4 + y4 * this.width;
//					bytes[index0] = (byte) tempData[y4][x4];
//				
//					this.transformedData[z - 1][y4][x4] = tempData[y4][x4];
//				}
//			}
//			ImagePlus imp = new ImagePlus("", new ByteProcessor(this.width, this.height, bytes));
//			imp.show();
//		}
//	
//		
//	}
	public static final byte[] float2Byte(float[] inData) {
		int j = 0;
		int length = inData.length;
		byte[] outData = new byte[length * 4];
		for (int i = 0; i < length; i++) {
			int data = Float.floatToIntBits(inData[i]);
			outData[j++] = (byte) (data >>> 24);
			outData[j++] = (byte) (data >>> 16);
			outData[j++] = (byte) (data >>> 8);
			outData[j++] = (byte) (data >>> 0);
		}
		return outData;
	}

	private void IWT(final double[] data) {
		final int dataLen = data.length;
		final int kernelLen = this.waveletFilter.rlf.length;
		final int mid = dataLen / 2;
		final int midKernel = kernelLen / 2;
		final int pad = dataLen * 20;
		final double[] approxUp = new double[dataLen];
		final double[] detailUp = new double[dataLen];
		final double[] temp = new double[dataLen];
		for (int i = 0; i < mid; ++i) {
			final int k = i * 2;
			approxUp[k] = data[i];
			approxUp[k + 1] = 0.0;
			detailUp[k] = data[mid + i];
			detailUp[k + 1] = 0.0;
		}
		for (int i = 0; i < dataLen; ++i) {
			double sumL = 0.0;
			double sumH = 0.0;
			for (int j = 0; j < kernelLen; ++j) {
				sumL += approxUp[(i - midKernel + j + pad) % dataLen] * this.waveletFilter.rlf[j];
				sumH += detailUp[(i - midKernel + j + pad) % dataLen] * this.waveletFilter.rhf[j];
			}
			temp[i] = sumL + sumH;
		}
		for (int i = 0; i < dataLen - 1; ++i) {
			data[i] = temp[i + 1];
		}
		data[dataLen - 1] = temp[0];
//		byte[] pixels = new byte[data.length];
//		for (int i = 0; i < data.length; i++)
//			pixels[i] = (byte) data[i];
//		IJ.log("pasaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
//		new ImagePlus("", new ByteProcessor(60, 50, pixels)).show();
//		IJ.log("pasaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa yjuuuuuuuuuuuuuuuuu");

	}

	private void IWT2D(final int z) {

		final double[][] tempData = new double[this.height][this.width];
		for (int y = 0; y < this.height; ++y) {
			for (int x = 0; x < this.width; ++x) {
				tempData[y][x] = this.transformedData[z - 1][y][x];
			}
		}
		for (int k = this.scale - 1; k >= 0; --k) {
			final int lev = 1 << k;
			final int levCols = this.height / lev;
			final int levRows = this.width / lev;
			final double[] col = new double[levRows];
			for (int y2 = 0; y2 < levCols; ++y2) {
				for (int x2 = 0; x2 < col.length; ++x2) {
					col[x2] = tempData[y2][x2];
				}
				this.IWT(col);

				for (int x2 = 0; x2 < col.length; ++x2) {
					tempData[y2][x2] = col[x2];

				}

			}
			final double[] row = new double[levCols];
			byte[] pixels = new byte[this.imageSize];
			for (int x3 = 0; x3 < levRows; ++x3) {
				for (int y3 = 0; y3 < row.length; ++y3) {
					row[y3] = tempData[y3][x3];
				}
				this.IWT(row);

				for (int y3 = 0; y3 < row.length; ++y3) {

					final int index = x3 + y3 * this.width;
					tempData[y3][x3] = row[y3];
					pixels[index] = (byte) row[y3];

				}

			}

		}
		for (int y4 = 0; y4 < this.height; ++y4) {
			for (int x4 = 0; x4 < this.width; ++x4) {
				this.imageData[z - 1][y4][x4] = tempData[y4][x4];

			}
		}

	}

	protected void fwdTransform(final int z) {
		this.FWT2D(z);
		this.prepareStretchWC(z);
	}

	protected void invTransform(final int z) {
		this.dataToWave(z);
		this.IWT2D(z);
	}

}