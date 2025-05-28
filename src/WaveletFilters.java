
import java.io.InputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import ij.IJ;

class WaveletFilters {
	protected double[] dlf;
	protected double[] dhf;
	protected double[] rlf;
	protected double[] rhf;

	protected void setFilter(final String filterName) {
		// final String content = this.getText(filterName);
		// if (!content.isEmpty()) {
		// final String[] lines = content.split("\\r?\\n|\\r");
		// final int filterSize = (lines.length - 3) / 4;
		this.dlf = new double[] { -0.12940952255092145, 0.22414386804185735, 0.836516303737469, 0.48296291314469025 };
		this.dhf = new double[] { -0.48296291314469025, 0.836516303737469, -0.22414386804185735, -0.12940952255092145 };
		this.rlf = new double[] { 0.48296291314469025, 0.836516303737469, 0.22414386804185735, -0.12940952255092145 };
		this.rhf = new double[] { -0.12940952255092145, -0.22414386804185735, 0.836516303737469, -0.48296291314469025 };
		int filterNumber = 0;
		int filterIndex = 0;

	}
}
