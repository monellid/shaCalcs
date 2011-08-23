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
	private double probExceed;
	private Site site;
	private EqkRupForecastAPI erf;
	private Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap;
	private List<Double> imlVals;
	private DisaggregationCalculator disCalc;

	@Before
	public void setUp() {
		latBinEdges = new double[] { -0.6, -0.4, -0.2, -0.0, 0.2, 0.4, 0.6 };
		lonBinEdges = new double[] { -0.6, -0.4, -0.2, -0.0, 0.2, 0.4, 0.6 };
		magBinEdges = new double[] { 5.0, 6.0, 7.0, 8.0, 9.0 };
		epsilonBinEdges = new double[] { -3.5, -2.5, -1.5, -0.5, +0.5, +1.5,
				+2.5, +3.5 };
		distanceBinEdges = new double[] { 60, 40, 20, 0 };
		probExceed = 0.1;
		site = new Site(new Location(0.0, 0.0));
		site.addParameter(new DoubleParameter(Vs30_Param.NAME, 760.0));
		erf = getTestERF();
		imrMap = new HashMap<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI>();
		imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, getTestIMR());
		imlVals = getTestImlVals();
		disCalc = new DisaggregationCalculator(latBinEdges, lonBinEdges,
				magBinEdges, epsilonBinEdges, distanceBinEdges);
		try {
			disCalc.disaggregate(probExceed, site, erf, imrMap, imlVals);
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
		double[] computedMagPMF = disCalc.getMagnitudePMF();

		Random rn = new Random(123456789);
		int n = 1000;

		// generate n stochastic event sets
		ArrayList<EqkRupture> ses = new ArrayList<EqkRupture>();
		for (int i = 0; i < n; i++) {
			ses.addAll(StochasticEventSetGenerator
					.getStochasticEventSetFromPoissonianERF(erf, rn));
		}

		// compute the ground motion value for which the disaggregation is
		// computed.
		double groundMotionValue;
		try {
			groundMotionValue = disCalc.computeGroundMotionValue(probExceed,
					site, erf, imrMap, imlVals);
		} catch (Exception e) {
			throw new RuntimeException(e.toString());
		}

		// for each magnitude bin, compute P(M=m | IML >x)
		double[] expectedMagPMF = new double[computedMagPMF.length];
		for (int i = 0; i < computedMagPMF.length; i++) {

			// loop over ruptures. If the rupture magnitude is in the considered
			// magnitude bin, compute a realization of the ground motion field
			// in the site of interest. Then check if the simulated ground
			// motion value is greater than then the ground motion value
			// corresponding to the probability of exceedence.
			double numRupture = 0.0;
			double totNumRupture = 0.0;
			List<Site> sites = new ArrayList<Site>();
			sites.add(site);
			for (EqkRupture rup : ses) {

				double mag = rup.getMag();

				if (mag >= magBinEdges[i] && mag < magBinEdges[i + 1]) {

					totNumRupture = totNumRupture + 1;

					ScalarIntensityMeasureRelationshipAPI attenRel = imrMap
							.get(rup.getTectRegType());
					GroundMotionFieldCalculator gmfCalc = new GroundMotionFieldCalculator(
							attenRel, rup, sites);
					Map<Site, Double> gmf = gmfCalc
							.getUncorrelatedGroundMotionField(rn);
					double gmfv = gmf.get(site);

					if (gmfv > groundMotionValue) {
						numRupture = numRupture + 1;
					}
				}

			}

			double expectedProbability = numRupture / totNumRupture;
			expectedMagPMF[i] = expectedProbability;
		}
		double totProb = 0;
		for(int i=0;i<expectedMagPMF.length;i++){
			totProb = totProb + expectedMagPMF[i];
		}
		for(int i=0;i<expectedMagPMF.length;i++){
			expectedMagPMF[i] = expectedMagPMF[i]/totProb;
		}
		
		for(int i=0;i<computedMagPMF.length;i++){
			System.out.println("mag: "+(magBinEdges[i]+magBinEdges[i+1])/2+", computed: "+computedMagPMF[i]+", expected: "+expectedMagPMF[i]);
		}

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

	private EqkRupForecastAPI getTestERF() {

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
