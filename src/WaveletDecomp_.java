import java.io.File;
import java.text.DecimalFormat;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Measurements;
import ij.plugin.PlugIn;
import ij.plugin.frame.SyncWindows;
import ij.process.ImageProcessor;

public class WaveletDecomp_ implements PlugIn, Measurements {
	ImagePlus[] imps;

	@Override
	public void run(String arg0) {
		File imagesDirectory = new File("/home/anaacayuela/Ana_pruebas_imageJ/javi_conesa/jonathan/images/canvas");
		File[] listOfFiles = imagesDirectory.listFiles();
		String[] imageTitles = new String[listOfFiles.length];
		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isFile())
				imageTitles[i] = listOfFiles[i].getName();
		}
		imps = new ImagePlus[imageTitles.length];
		for (int i = 0; i < imageTitles.length; i++) {
			imps[i] = IJ.openImage(imagesDirectory.getAbsolutePath() + File.separator + imageTitles[i]);
			// ImagePlus[] slices = stack2images(imps[i]);
			// for (int j = 0; j < slices.length; j++) {

			ImageData imageData = new ImageData(imps[i]);
			String title = imageData.imageOrig.getTitle();
			String newTitle = "WT-" + title;
			imageData.imageWave.setTitle(newTitle);
			newTitle = "WT-NoStretch-" + title;
//			imageData.imageWaveNoStretch.setTitle(newTitle);
			title = imageData.imageOrig.getTitle();
			newTitle = "Filtered-" + title;
			imageData.imageModif.setTitle(newTitle);
			// imageData.imageModif.show();
			int maxScale = imageData.changeScale(1);
			imageData.changeScale(1);
			imageData.setScale(8);
			imageData.setWaveletFilter("Daubechies 2");

			// transform
			for (int z = 1; z <= imageData.getNSlices(); ++z) {
				imageData.bitmapToData(z);
				imageData.fwdTransform(z);
				imageData.invTransform(z);
				imageData.dataToBitmap(z);
			}

			imageData.imageWave.show();
			imageData.imageWave.updateAndDraw();
			// }
			// imageData.imageModif.show();
			imageData.imageModif.updateAndDraw();

			// if (Wavelet_Denoise.showNoStretchWT) {
			imageData.imageWaveNoStretch.show();
			imageData.imageWaveNoStretch.updateAndDraw();

			// }
			synchronizeImages();
			// }
		}

	}

	public static ImagePlus[] stack2images(ImagePlus imp) {
		String sLabel = imp.getTitle();
		String sImLabel = "";
		ImageStack stack = imp.getStack();

		int sz = stack.getSize();
		int currentSlice = imp.getCurrentSlice(); // to reset ***

		DecimalFormat df = new DecimalFormat("0000"); // for title
		ImagePlus[] arrayOfImages = new ImagePlus[imp.getStack().getSize()];
		for (int n = 1; n <= sz; ++n) {
			imp.setSlice(n); // activate next slice ***

			// Get current image processor from stack. What ever is
			// used here should do a COPY pixels from old processor to
			// new. For instance, ImageProcessor.crop() returns copy.
			ImageProcessor ip = imp.getProcessor(); // ***
			ImageProcessor newip = ip.createProcessor(ip.getWidth(), ip.getHeight());
			newip.setPixels(ip.getPixelsCopy());

			// Create a suitable label, using the slice label if possible
			sImLabel = imp.getStack().getSliceLabel(n);
			if (sImLabel == null || sImLabel.length() < 1) {
				sImLabel = "slice" + df.format(n) + "_" + sLabel;
			}
			// Create new image corresponding to this slice.
			ImagePlus im = new ImagePlus(sImLabel, newip);
			im.setCalibration(imp.getCalibration());
			arrayOfImages[n - 1] = im;

			// Show this image.
			// imp.show();
		}
		// Reset original stack state.
		imp.setSlice(currentSlice); // ***
		if (imp.isProcessor()) {
			ImageProcessor ip = imp.getProcessor();
			ip.setPixels(ip.getPixels()); // ***
		}
		imp.setSlice(currentSlice);
		return arrayOfImages;
	}

	private enum Mode {
		Nothing("Nothing", 0), Suppress("Suppress", 1), SoftDenoise("SoftDenoise", 2), HardDenoise("HardDenoise", 3),
		SuppressApprox("SuppressApprox", 4), SuppressDetail("SuppressDetail", 5);

		private Mode(final String name, final int ordinal) {
		}
	}

	private void synchronizeImages() {
		final SyncWindows syncW = new SyncWindows();
		syncW.run((String) null);
	}

}
