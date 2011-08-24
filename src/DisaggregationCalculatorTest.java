import static org.junit.Assert.*;

import java.rmi.RemoteException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.param.DoubleParameter;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.griddedForecast.MagFreqDistsForFocalMechs;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.GEM1ERF;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.SourceData.GEMAreaSourceData;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.SourceData.GEMSourceData;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.attenRelImpl.BA_2008_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.OtherParams.ComponentParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

public class DisaggregationCalculatorTest {

	private double[] latBinEdges;
	private double[] lonBinEdges;
	private double[] magBinEdges;
	private double[] epsilonBinEdges;
	private double[] distanceBinEdges;
	private double probExceed = 0.1;
	private Site site;
	private GEM1ERF erf;
	private Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap;
	private List<Double> imlVals;
	private DisaggregationCalculator disCalc;

	private static Random rn = new Random(123456789);
	private static int n = 2000;
	private ArrayList<EqkRupture> ses;
	private double groundMotionValue;

	@Before
	public void setUp() {
		latBinEdges = new double[] { -0.6, -0.4, -0.2, -0.0, 0.2, 0.4, 0.6 };
		lonBinEdges = new double[] { -0.6, -0.4, -0.2, -0.0, 0.2, 0.4, 0.6 };
		magBinEdges = new double[] { 5.0, 6.0, 7.0, 8.0, 9.0 };
		epsilonBinEdges = new double[] { -6.5, -5.5, -4.5, -3.5, -2.5, -1.5,
				-0.5, +0.5, +1.5, +2.5, +3.5, +4.5, +5.5, +6.5 };
		distanceBinEdges = new double[] { 60, 40, 20, 0 };
		site = new Site(new Location(0.0, 0.0));
		site.addParameter(new DoubleParameter(Vs30_Param.NAME, 760.0));
		erf = getTestERF();
		imrMap = new HashMap<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI>();
		imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, getTestIMR());
		imlVals = getTestImlVals();

		// compute disaggregation matrix
		disCalc = new DisaggregationCalculator(latBinEdges, lonBinEdges,
				magBinEdges, epsilonBinEdges, distanceBinEdges);
		try {
			disCalc.disaggregate(probExceed, site, erf, imrMap, imlVals);
		} catch (Exception e) {
			throw new RuntimeException(e.toString());
		}
		ses = new ArrayList<EqkRupture>();

		// generate n stochastic event sets
		for (int i = 0; i < n; i++) {
			ses.addAll(StochasticEventSetGenerator
					.getStochasticEventSetFromPoissonianERF(erf, rn));
		}

		// compute the ground motion value for which the disaggregation is
		// computed.
		try {
			groundMotionValue = disCalc.computeGroundMotionValue(probExceed,
					site, erf, imrMap, imlVals);
		} catch (Exception e) {
			throw new RuntimeException(e.toString());
		}
	}

	@After
	public void tearDown() {
		disCalc = null;
	}

	/**
	 * The magnitude PMF is checked against the magnitude PMF computed using
	 * Monte Carlo approach.
	 * 
	 */
	@Test
	public void getMagnitudePMFTest() {

		// computed mag PMF
		double[] computedMagPMF = disCalc.getMagnitudePMF();

		double[] expectedMagPMF = compute1DPMFThroughMonteCarlo("magnitude",
				magBinEdges);

		for (int i = 0; i < computedMagPMF.length; i++) {
			System.out.println("mag: " + (magBinEdges[i] + magBinEdges[i + 1])
					/ 2 + ", computed: " + computedMagPMF[i] + ", expected: "
					+ expectedMagPMF[i]);
			assertEquals(computedMagPMF[i], expectedMagPMF[i], 0.1);
		}

	}

	// for each magnitude bin, compute:
	// P(M=m | IML >x) = lambda(IML > x, M =m) / lambda (IML > x)
	// loop over ruptures. For each rupture, compute ground motion field
	// realization in the site of interest. Then check if the simulated
	// ground motion field gives a value greater than the ground motion
	// value corresponding to the probability of interest.
	// lambda(IML > x, M =m) is given by the number of ruptures with
	// magnitude inside the considered magnitude bin, and that produce a
	// ground motion value greater then the ground motion value of interest.
	// lambda (IML > x) is the number of ruptures which produce a ground
	// motion value greater than the one of interest (indipendently of the
	// magnitude)
	private double[] compute1DPMFThroughMonteCarlo(String paramName,
			double[] paramBinEdges) {
		double[] pmf = new double[paramBinEdges.length - 1];
		List<Site> sites = new ArrayList<Site>();
		sites.add(site);
		for (int i = 0; i < pmf.length; i++) {

			double numRuptureWithGivenM = 0.0;
			double numRupture = 0.0;

			for (EqkRupture rup : ses) {

				// check if the rupture is inside the ranges considered
				// find closest point in the rupture area
				Location closestLoc = disCalc.getClosestLocation(site,
						rup.getRuptureSurface());

				ScalarIntensityMeasureRelationshipAPI imr = imrMap.get(rup
						.getTectRegType());
				imr.setSite(site);
				imr.setEqkRupture(rup);
				imr.setIntensityMeasureLevel(groundMotionValue);

				double lat = closestLoc.getLatitude();
				double lon = closestLoc.getLongitude();
				double magnitude = rup.getMag();
				double epsilon = imr.getEpsilon();
				String tectonicRegionType = rup.getTectRegType().toString();

				// if one of the parameters (latitude, longitude, magnitude,
				// epsilon) is outside of
				// the considered range do not include in the conditional
				// probability calculation, that is skip the rupture
				if (lat < latBinEdges[0]
						|| lat >= latBinEdges[latBinEdges.length - 1]) {
					continue;
				}
				if (lon < lonBinEdges[0]
						|| lon >= lonBinEdges[lonBinEdges.length - 1]) {
					continue;
				}
				if (magnitude < magBinEdges[0]
						|| magnitude >= magBinEdges[magBinEdges.length - 1]) {
					continue;
				}
				if (epsilon < epsilonBinEdges[0]
						|| epsilon >= epsilonBinEdges[epsilonBinEdges.length - 1]) {
					continue;
				}

				// compute ground motion field
				ScalarIntensityMeasureRelationshipAPI attenRel = imrMap.get(rup
						.getTectRegType());
				GroundMotionFieldCalculator gmfCalc = new GroundMotionFieldCalculator(
						attenRel, rup, sites);
				Map<Site, Double> gmf = gmfCalc
						.getUncorrelatedGroundMotionField(rn);
				double gmfv = gmf.get(site);

				// check that the value is greater than the value of interest
				if (gmfv > groundMotionValue) {
					numRupture = numRupture + 1;

					double param;
					if (paramName.equalsIgnoreCase("magnitude")) {
						param = magnitude;
					} else {
						throw new RuntimeException(
								"parameter name not recognized!!");
					}
					if (param >= paramBinEdges[i]
							&& param < paramBinEdges[i + 1]) {
						numRuptureWithGivenM = numRuptureWithGivenM + 1;
					}
				}
			}

			double expectedProbability = numRuptureWithGivenM / numRupture;
			pmf[i] = expectedProbability;
		}
		return pmf;
	}

	private List<Double> getTestImlVals() {

		List<Double> imlVals = new ArrayList<Double>();
		imlVals.add(Math.log(0.005));
		imlVals.add(Math.log(0.007));
		imlVals.add(Math.log(0.0098));
		imlVals.add(Math.log(0.0137));
		imlVals.add(Math.log(0.0192));
		imlVals.add(Math.log(0.0269));
		imlVals.add(Math.log(0.0376));
		imlVals.add(Math.log(0.0527));
		imlVals.add(Math.log(0.0738));
		imlVals.add(Math.log(0.103));
		imlVals.add(Math.log(0.145));
		imlVals.add(Math.log(0.203));
		imlVals.add(Math.log(0.284));
		imlVals.add(Math.log(0.397));
		imlVals.add(Math.log(0.556));
		imlVals.add(Math.log(0.778));
		imlVals.add(Math.log(1.09));
		imlVals.add(Math.log(1.52));
		imlVals.add(Math.log(2.13));
		return imlVals;
	}

	private ScalarIntensityMeasureRelationshipAPI getTestIMR() {
		ScalarIntensityMeasureRelationshipAPI imr = new BA_2008_AttenRel(null);
		imr.setIntensityMeasure(PGA_Param.NAME);
		imr.getParameter(StdDevTypeParam.NAME).setValue(
				StdDevTypeParam.STD_DEV_TYPE_TOTAL);
		imr.getParameter(SigmaTruncTypeParam.NAME).setValue(
				SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_2SIDED);
		imr.getParameter(SigmaTruncLevelParam.NAME).setValue(3.0);
		imr.getParameter(ComponentParam.NAME).setValue(
				ComponentParam.COMPONENT_GMRotI50);
		return imr;
	}

	private GEM1ERF getTestERF() {

		ArrayList<GEMSourceData> srcList = new ArrayList<GEMSourceData>();
		srcList.add(getTestSourceData());
		double timeSpan = 50.0;
		return GEM1ERF.getGEM1ERF(srcList, timeSpan);
	}

	private GEMSourceData getTestSourceData() {

		String id = "src1";
		String name = "testSource";
		TectonicRegionType tectReg = TectonicRegionType.ACTIVE_SHALLOW;
		LocationList border = new LocationList();
		border.add(new Location(-0.5, -0.5));
		border.add(new Location(0.5, -0.5));
		border.add(new Location(0.5, 0.5));
		border.add(new Location(-0.5, 0.5));
		Region reg = new Region(border, BorderType.GREAT_CIRCLE);
		double bValue = 1.0;
		double totCumRate = 0.2;
		double min = 5.05;
		double max = 8.95;
		int num = 41;
		GutenbergRichterMagFreqDist magDist = new GutenbergRichterMagFreqDist(
				bValue, totCumRate, min, max, num);
		double strike = 0.0;
		double dip = 90.0;
		double rake = 0.0;
		FocalMechanism focalMechanism = new FocalMechanism(strike, dip, rake);
		MagFreqDistsForFocalMechs magfreqDistFocMech = new MagFreqDistsForFocalMechs(
				magDist, focalMechanism);
		ArbitrarilyDiscretizedFunc aveRupTopVsMag = new ArbitrarilyDiscretizedFunc();
		double magThreshold = 6.5;
		double topOfRuptureDepth = 0.0;
		aveRupTopVsMag.set(magThreshold, topOfRuptureDepth);
		double aveHypoDepth = 5.0;
		GEMSourceData srcData = new GEMAreaSourceData(id, name, tectReg, reg,
				magfreqDistFocMech, aveRupTopVsMag, aveHypoDepth);
		return srcData;
	}
}
