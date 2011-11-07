import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.TectonicRegionType;

public class UniformHazardSpectrumCalculatorTest {

	private double[] periods;
	private EqkRupForecastAPI erf;
	private Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap;
	private List<Double> imlVals;
	private double maxDistance;

	@Before
	public void setUp() {
		periods = CalculatorsTestHelper.getTestIMRPeriods();
		erf = CalculatorsTestHelper.getTestERF();
		imrMap = CalculatorsTestHelper.getTestImrMap();
		imlVals = CalculatorsTestHelper.getTestImlVals();
		maxDistance = 200.0;

	}

	@After
	public void tearDown() {
		periods = null;
		erf = null;
		imrMap = null;
		imlVals = null;
	}

	/**
	 * Checks UHS calculation on supported SA periods. The test is done by
	 * comparing the calculated UHS, with one obtained by calculating the hazard
	 * curve for each period, and then extracting from each hazard curve the
	 * ground motion value corresponding to the period of interest.
	 */
	@Test
	public void UHSOnSupportedSAPeriodsTest() {

		double[] poE = new double[]{0.1, 0.02};

		UniformHazardSpectrumCalculator uhsCalc = new UniformHazardSpectrumCalculator(
				periods, erf, imrMap, imlVals, maxDistance);

		List<double[]> calculatedUHS = uhsCalc.getUHS(poE,
				CalculatorsTestHelper.getTestSite());

		List<double[]> expectedUHS = computeExpectedUHSForSupportedSA(poE);

		for(int ip=0;ip<poE.length;ip++){
			for (int i = 0; i < periods.length; i++) {
				assertEquals(calculatedUHS.get(ip)[i], expectedUHS.get(ip)[i], 1e-5);
			}
		}
	}

	/**
	 * Checks UHS calculation on PGA and supported SA periods.
	 */
	@Test
	public void UHSOnPGAAndSupportedSAPeriodsTest() {

		double[] poE = new double[]{0.1,0.02};

		double[] periods = new double[this.periods.length + 1];
		periods[0] = 0.0;
		for (int i = 1; i < periods.length; i++) {
			periods[i] = this.periods[i - 1];
		}

		UniformHazardSpectrumCalculator uhsCalc = new UniformHazardSpectrumCalculator(
				periods, erf, imrMap, imlVals, maxDistance);

		List<double[]> calculatedUHS = uhsCalc.getUHS(poE,
				CalculatorsTestHelper.getTestSite());

		List<double[]> expectedUHS = computeExpectedUHSForPGAAndSupportedSA(poE);

		for(int ip=0;ip<poE.length;ip++){
			System.out.println("period: "+poE[ip]);
			for (int i = 0; i < periods.length; i++) {
				System.out.println("calculated UHS: "+Math.exp(calculatedUHS.get(ip)[i])+", expected: "+Math.exp(expectedUHS.get(ip)[i]));
				assertEquals(calculatedUHS.get(ip)[i], expectedUHS.get(ip)[i], 1e-5);
			}
			System.out.println("\n");
		}

	}

	/**
	 * Checks UHS calculation on not-supported SA periods, requiring therefore
	 * interpolation. The reference UHS (expected) is calculated by calculating
	 * a UHS on the supported periods that are closer to the requested ones, and
	 * then by interpolating this UHS on the requested periods (using linear
	 * interpolation).
	 */
	@Test
	public void UHSOnInterpolatedSAPeriodsTest() {

		double[] poE = new double[]{0.1, 0.02};

		// these periods are not supported by the BA2008 gmpe used for the test
		// and therefore require interpolation
		double[] periods = new double[] { 0.025, 0.45, 2.5 };

		UniformHazardSpectrumCalculator uhsCalc = new UniformHazardSpectrumCalculator(
				periods, erf, imrMap, imlVals, maxDistance);

		List<double[]> calculatedUHS = uhsCalc.getUHS(poE,
				CalculatorsTestHelper.getTestSite());

		List<double[]> expectedUHS = computeExpectedUHSForInterpolatedSA(poE);

		for(int ip=0;ip<poE.length;ip++){
			for (int i = 0; i < periods.length; i++) {
				assertEquals(calculatedUHS.get(ip)[i], expectedUHS.get(ip)[i], 1e-1);
			}	
		}

	}

	private List<double[]> computeExpectedUHSForSupportedSA(double[] poE) {
		List<double[]> expectedUHSSet = new ArrayList<double[]>();
		for(int i=0;i<poE.length;i++){
			double[] expectedUHS = new double[periods.length];
			int periodIndex = 0;
			for (double period : periods) {
				setImrPeriod(period);
				expectedUHS[periodIndex] = CalculatorsUtils
						.computeGroundMotionValue(poE[i],
								CalculatorsTestHelper.getTestSite(), erf, imrMap,
								imlVals);
				periodIndex = periodIndex + 1;
			}
			expectedUHSSet.add(expectedUHS);
		}

		return expectedUHSSet;
	}

	private List<double[]> computeExpectedUHSForPGAAndSupportedSA(double[] poE) {

		List<double[]> expectedUHSSet = new ArrayList<double[]>();
		
		for(int i=0;i<poE.length;i++){
			
			double[] expectedUHS = new double[periods.length + 1];
			
			// compute UHS for PGA
			for (ScalarIntensityMeasureRelationshipAPI imr : imrMap.values()) {
				imr.setIntensityMeasure(PGA_Param.NAME);
			}
			expectedUHS[0] = CalculatorsUtils.computeGroundMotionValue(poE[i],
					CalculatorsTestHelper.getTestSite(), erf, imrMap, imlVals);

			// compute UHS for SA periods
			int periodIndex = 1;
			for (double period : periods) {
				setImrPeriod(period);
				expectedUHS[periodIndex] = CalculatorsUtils
						.computeGroundMotionValue(poE[i],
								CalculatorsTestHelper.getTestSite(), erf, imrMap,
								imlVals);
				periodIndex = periodIndex + 1;
			}
			
			expectedUHSSet.add(expectedUHS);
		}
		
		return expectedUHSSet;
	}

	private List<double[]> computeExpectedUHSForInterpolatedSA(double[] poE) {

		List<double[]> expectedUHSSet = new ArrayList<double[]>();
		
		for(int i=0;i<poE.length;i++){
			
			double[] expectedUHS = new double[periods.length];

			double period = 0.025;
			double period1 = 0.02;
			double period2 = 0.03;
			double gmv = getInterpolatedUHSValue(poE[i], period, period1, period2);
			expectedUHS[0] = gmv;

			period = 0.45;
			period1 = 0.4;
			period2 = 0.5;
			gmv = getInterpolatedUHSValue(poE[i], period, period1, period2);
			expectedUHS[1] = gmv;

			period = 2.5;
			period1 = 2.0;
			period2 = 3.0;
			gmv = getInterpolatedUHSValue(poE[i], period, period1, period2);
			expectedUHS[2] = gmv;	
			
			expectedUHSSet.add(expectedUHS);
		}

		return expectedUHSSet;
	}

	private double getInterpolatedUHSValue(double poE, double period,
			double period1, double period2) {

		setImrPeriod(period1);
		double gmv1 = CalculatorsUtils.computeGroundMotionValue(poE,
				CalculatorsTestHelper.getTestSite(), erf, imrMap, imlVals);
		setImrPeriod(period2);
		double gmv2 = CalculatorsUtils.computeGroundMotionValue(poE,
				CalculatorsTestHelper.getTestSite(), erf, imrMap, imlVals);

		double A = (period2 - period) / (period2 - period1);
		double B = 1 - A;
		double gmv = A * gmv1 + B * gmv2;
		return gmv;
	}

	private void setImrPeriod(double period) {
		for (ScalarIntensityMeasureRelationshipAPI imr : imrMap.values()) {
			imr.setIntensityMeasure(SA_Param.NAME);
			imr.getParameter(PeriodParam.NAME).setValue(period);
		}
	}

}
