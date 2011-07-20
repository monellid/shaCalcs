import java.rmi.RemoteException;
import java.util.HashMap;
import java.util.Map;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.StringParameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.earthquake.rupForecastImpl.PoissonianAreaSourceTestHelper;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.util.TectonicRegionType;

public class DisaggregationCalculatorTest {

	private static EqkRupForecastAPI erf;
	private static ScalarIntensityMeasureRelationshipAPI imr;

	@Before
	public void setUp() {
		// area source: min mag = 5.0, max mag = 6.5.
		erf = PoissonianAreaSourceTestHelper
				.getPeerTestSet1Case10AreaSourceErf();
		// Sadigh et Al. 1997, PGA, no standard devitation
		imr = PoissonianAreaSourceTestHelper.getPeerTestSet1Case10GMPE();
		// set standard deviation to total
		imr.getParameter(StdDevTypeParam.NAME).setValue(
				StdDevTypeParam.STD_DEV_TYPE_TOTAL);
	}

	@After
	public void tearDown() {
		erf = null;
		imr = null;
	}

	@Test
	public void disaggregationTest() throws RemoteException {

		Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap = new HashMap<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI>();
		imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, imr);

		// compute hazard curve for middle point
		StringParameter sadighSiteType = new StringParameter("Sadigh Site Type");
		sadighSiteType.setValue("Rock");
		Site site = new Site(new Location(38.000, -122.000));
		site.addParameter(sadighSiteType);

		// hazard curve
		ArbitrarilyDiscretizedFunc hazCurve = new ArbitrarilyDiscretizedFunc();
		hazCurve.set(Math.log(0.001), 1.0);
		hazCurve.set(Math.log(0.01), 1.0);
		hazCurve.set(Math.log(0.05), 1.0);
		hazCurve.set(Math.log(0.1), 1.0);
		hazCurve.set(Math.log(0.15), 1.0);
		hazCurve.set(Math.log(0.2), 1.0);
		hazCurve.set(Math.log(0.25), 1.0);
		hazCurve.set(Math.log(0.3), 1.0);
		hazCurve.set(Math.log(0.35), 1.0);
		hazCurve.set(Math.log(0.4), 1.0);

		// compute hazard curve
		HazardCurveCalculator hcc = new HazardCurveCalculator();
		hcc.getHazardCurve(hazCurve, site, imrMap, erf);

		// compute ground motion corresponding to 1e-4 probability of exceedence
		double probLevel = 1e-4;
		System.out.println("hazard curve: " + hazCurve);
		double iml = hazCurve.getFirstInterpolatedX(probLevel);

		// compute disaggregation
		double[] mag = new double[] { 5.25,5.75,6.25,6.75};
		double[] dist = new double[] { 5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0,
				75.0, 85.0, 95.0, 105.0 };
		double[] epsilon = new double[] { -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0 };
		DisaggregationCalculator disaggCalc = new DisaggregationCalculator(mag,
				dist, epsilon);
		double[][][] pmf = disaggCalc.getMagDistEpsilonDisaggregation(iml,
				site, erf, imrMap);

		for (int i = 0; i < epsilon.length; i++) {
			System.out.println("Epsilon: " + epsilon[i]);
			System.out.print("      ");
			for (int k = 0; k < dist.length; k++) {
				System.out.print(String.format("%2.2f    ", dist[k]));
			}
			System.out.println();
			for (int j = 0; j < mag.length; j++) {
				System.out.print(String.format("%2.2f ", mag[j]));
				for (int k = 0; k < dist.length; k++) {
					System.out.print(String.format("%.1e ", pmf[j][k][i]));
					// compare against 
				}
				System.out.println();
			}
			System.out.println();
		}

		// sum all the PMFs contributions
		double totSum = 0.0;
		for (int i = 0; i < mag.length; i++) {
			for (int j = 0; j < dist.length; j++) {
				for (int k = 0; k < epsilon.length; k++) {
					totSum = totSum + pmf[i][j][k];
				}
			}
		}
		System.out.println("Sum of all contributions: " + totSum);

	}

}