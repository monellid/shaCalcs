import static org.junit.Assert.*;

import java.rmi.RemoteException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.BlockFieldMatrix;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
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
	private static int n = 3000;
	private DisaggregationCalculatorMC disCalcMC;

	@Before
	public void setUp() {
		// initialize input
		latBinEdges = new double[] {-0.6,-0.3,-0.1,0.1,0.3,0.6};//{ -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6 };
		lonBinEdges = new double[] {-0.6,-0.3,-0.1,0.1,0.3,0.6};//{ -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6 };
		magBinEdges = new double[] { 5.0, 6.0, 7.0, 8.0, 9.0 };
		epsilonBinEdges = new double[] {-3.5, -2.5, -1.5,
				-0.5, +0.5, +1.5, +2.5, +3.5};
		distanceBinEdges = new double[] { 0, 20, 40, 60 };
		site = new Site(new Location(0.0, 0.0));
		site.addParameter(new DoubleParameter(Vs30_Param.NAME, 760.0));
		erf = getTestERF();
		imrMap = new HashMap<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI>();
		imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, getTestIMR());
		imlVals = getTestImlVals();

		// compute disaggregation matrix following classical approach
		disCalc = new DisaggregationCalculator(latBinEdges, lonBinEdges,
				magBinEdges, epsilonBinEdges, distanceBinEdges);
		disCalc.disaggregate(probExceed, site, erf, imrMap, imlVals);

		// compute disaggregation calculator with MC approach
		disCalcMC = new DisaggregationCalculatorMC(latBinEdges, lonBinEdges,
				magBinEdges, epsilonBinEdges, distanceBinEdges);
		disCalcMC.disaggregate(probExceed, site, erf, imrMap, imlVals, rn, n);
	}

	@After
	public void tearDown() {
		disCalc = null;
		disCalcMC = null;
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

		double[] expectedMagPMF = disCalcMC.getMagnitudePMF();

		for (int i = 0; i < computedMagPMF.length; i++) {
			System.out.println("mag: " + (magBinEdges[i] + magBinEdges[i + 1])
					/ 2 + ", computed: " + computedMagPMF[i] + ", expected: "
					+ expectedMagPMF[i]);
			assertEquals(computedMagPMF[i], expectedMagPMF[i], 0.1);
		}

	}

	/**
	 * The distance PMF is checked against the distance PMF computed using Monte
	 * Carlo approach.
	 */
	@Test
	public void getDistancePMFTest() {

		double[] computedDistancePMF = disCalc.getDistancePMF();

		double[] expectedDistancePMF = disCalcMC.getDistancePMF();

		for (int i = 0; i < computedDistancePMF.length; i++) {
			System.out.println("distance (km): "
					+ (distanceBinEdges[i] + distanceBinEdges[i + 1]) / 2
					+ ", computed: " + computedDistancePMF[i] + ", expected: "
					+ expectedDistancePMF[i]);
			assertEquals(computedDistancePMF[i], expectedDistancePMF[i], 0.1);
		}
	}

	/**
	 * The tectonic region type PMF is checked against the tectonic region type
	 * PMF computed using Monte Carlo approach.
	 */
	@Test
	public void getTectonicRegionTypePMFTest() {

		double[] computedTectonicRegionTypePMF = disCalc
				.getTectonicRegionTypePMF();

		double[] expectedTectonicRegionTypePMF = disCalcMC
				.getTectonicRegionTypePMF();
		
		for (int i = 0; i < computedTectonicRegionTypePMF.length; i++) {
			System.out.println("tectonic region type: "+(i+1)
					+ ", computed: " + computedTectonicRegionTypePMF[i] + ", expected: "
					+ expectedTectonicRegionTypePMF[i]);
			assertEquals(computedTectonicRegionTypePMF[i], expectedTectonicRegionTypePMF[i], 0.1);
		}
	}
	
	@Test
	public void getMagnitudeTectonicRegionTypePMFTest(){
		double[][] computedMagTrtPMF = disCalc.getMagnitudeTectonicRegionTypePMF();
		
		double[][] expectedMagTrtPMF = disCalcMC.getMagnitudeTectonicRegionTypePMF();
		
		for(int i=0;i<computedMagTrtPMF.length;i++){
			for(int j=0;j<computedMagTrtPMF[i].length;j++){
				assertEquals(computedMagTrtPMF[i][j], expectedMagTrtPMF[i][j], 0.1);
			}
		}
		
	}
	
	
	@Test
	public void getMagnitudeDistancePMFTest(){
		double[][] computedMagDistPMF = disCalc.getMagnitudeDistancePMF();
		
		double[][] expectedMagDistPMF = disCalcMC.getMagnitudeDistancePMF();
		
		for(int i=0;i<computedMagDistPMF.length;i++){
			for(int j=0;j<computedMagDistPMF[i].length;j++){
				assertEquals(computedMagDistPMF[i][j], expectedMagDistPMF[i][j], 0.1);
			}
		}
		
		for(int j=0;j<computedMagDistPMF[0].length;j++){
			System.out.print((distanceBinEdges[j]+distanceBinEdges[j+1])/2+" ");
		}
		System.out.println();
		for(int i=0;i<computedMagDistPMF.length;i++){
			System.out.print((magBinEdges[i]+magBinEdges[i+1])/2+" ");
			for(int j=0;j<computedMagDistPMF[i].length;j++){
				System.out.print(computedMagDistPMF[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	@Test
	public void getMagnitudeDistanceEpsilonPMFTest(){
		double[][][] computedMagDistEpsilonPMF = disCalc.getMagnitudeDistanceEpsilonPMF();
		
		double[][][] expectedMagDistEpsilonPMF = disCalcMC.getMagnitudeDistanceEpsilonPMF();
		
		for(int i=0;i<computedMagDistEpsilonPMF.length;i++){
			for(int j=0;j<computedMagDistEpsilonPMF[i].length;j++){
				for(int k=0;k<computedMagDistEpsilonPMF[i][j].length;k++){
					assertEquals(computedMagDistEpsilonPMF[i][j][k], expectedMagDistEpsilonPMF[i][j][k], 0.1);
				}
			}
		}
		for(int k=0;k<computedMagDistEpsilonPMF[0][0].length;k++){
			System.out.println("Epsilon: "+(epsilonBinEdges[k]+epsilonBinEdges[k+1])/2);
			for(int j=0;j<computedMagDistEpsilonPMF[0].length;j++){
				System.out.print((distanceBinEdges[j]+distanceBinEdges[j+1])/2+" ");
			}
			System.out.println();
			for(int i=0;i<computedMagDistEpsilonPMF.length;i++){
				System.out.print((magBinEdges[i]+magBinEdges[i+1])/2+" ");
				for(int j=0;j<computedMagDistEpsilonPMF[i].length;j++){
					System.out.print(computedMagDistEpsilonPMF[i][j][k]+" ");
				}
				System.out.println();
			}
		}
	}
	
	@Test
	public void getLatituteLongitudePMFTest(){
		double[][] computedLatLonPMF = disCalc.getLatitudeLongitudePMF();
		double[][] expectedLatLonPMF = disCalcMC.getLatitudeLongitudePMF();
		
		for(int i=0;i<computedLatLonPMF.length;i++){
			for(int j=0;j<computedLatLonPMF[i].length;j++){
				assertEquals(computedLatLonPMF[i][j], expectedLatLonPMF[i][j], 0.1);
			}
		}
		
		for(int j=0;j<computedLatLonPMF[0].length;j++){
			System.out.print((latBinEdges[j]+latBinEdges[j+1])/2+" ");
		}
		System.out.println();
		for(int i=0;i<computedLatLonPMF.length;i++){
			System.out.print((lonBinEdges[i]+lonBinEdges[i+1])/2+" ");
			for(int j=0;j<computedLatLonPMF[i].length;j++){
				System.out.print(computedLatLonPMF[i][j]+" ");
			}
			System.out.println();
		}
	}

	@Test
	public void getLatitudeLongitudeMagnitudePMFTest(){
		double[][][] computedLatLonMagPMF = disCalc.getLatitudeLongitudeMagnitudePMF();
		double[][][] expectedLatLonMagPMF = disCalcMC.getLatitudeLongitudeMagnitudePMF();
		
		for(int i=0;i<computedLatLonMagPMF.length;i++){
			for(int j=0;j<computedLatLonMagPMF[i].length;j++){
				for(int k=0;k<computedLatLonMagPMF[i][j].length;k++){
					assertEquals(computedLatLonMagPMF[i][j][k], expectedLatLonMagPMF[i][j][k], 0.1);
				}
			}
		}
		for(int k=0;k<computedLatLonMagPMF[0][0].length;k++){
			System.out.println("Magnitude: "+(magBinEdges[k]+magBinEdges[k+1])/2);
			for(int j=0;j<computedLatLonMagPMF[0].length;j++){
				System.out.print((lonBinEdges[j]+lonBinEdges[j+1])/2+" ");
			}
			System.out.println();
			for(int i=0;i<computedLatLonMagPMF.length;i++){
				System.out.print((latBinEdges[i]+latBinEdges[i+1])/2+" ");
				for(int j=0;j<computedLatLonMagPMF[i].length;j++){
					System.out.print(computedLatLonMagPMF[i][j][k]+" ");
				}
				System.out.println();
			}
		}
	}
	
	@Test
	public void getLatitudeLongitudeEpsilonPMFTest(){
		double[][][] computedLatLonEpsilonPMF = disCalc.getLatitudeLongitudeEpsilonPMF();
		double[][][] expectedLatLonEpsilonPMF = disCalcMC.getLatitudeLongitudeEpsilonPMF();
		
		for(int i=0;i<computedLatLonEpsilonPMF.length;i++){
			for(int j=0;j<computedLatLonEpsilonPMF[i].length;j++){
				for(int k=0;k<computedLatLonEpsilonPMF[i][j].length;k++){
					assertEquals(computedLatLonEpsilonPMF[i][j][k], expectedLatLonEpsilonPMF[i][j][k], 0.1);
				}
			}
		}
		for(int k=0;k<computedLatLonEpsilonPMF[0][0].length;k++){
			System.out.println("Epsilon: "+(epsilonBinEdges[k]+epsilonBinEdges[k+1])/2);
			for(int j=0;j<computedLatLonEpsilonPMF[0].length;j++){
				System.out.print((lonBinEdges[j]+lonBinEdges[j+1])/2+" ");
			}
			System.out.println();
			for(int i=0;i<computedLatLonEpsilonPMF.length;i++){
				System.out.print((latBinEdges[i]+latBinEdges[i+1])/2+" ");
				for(int j=0;j<computedLatLonEpsilonPMF[i].length;j++){
					System.out.print(computedLatLonEpsilonPMF[i][j][k]+" ");
				}
				System.out.println();
			}
		}
	}
	
	@Test
	public void getLatitudeLongitudeMagnitudeEpsilonPMFTest(){
		double[][][][] computedLatLonMagEpsilonPMF = disCalc.getLatitudeLongitudeMagnitudeEpsilonPMF();
		double[][][][] expectedLatLonMagEpsilonPMF = disCalcMC.getLatitudeLongitudeMagnitudeEpsilonPMF();
		
		for(int i=0;i<computedLatLonMagEpsilonPMF.length;i++){
			for(int j=0;j<computedLatLonMagEpsilonPMF[i].length;j++){
				for(int k=0;k<computedLatLonMagEpsilonPMF[i][j].length;k++){
					for(int l=0;l<computedLatLonMagEpsilonPMF[i][j][k].length;l++){
						assertEquals(computedLatLonMagEpsilonPMF[i][j][k][l], expectedLatLonMagEpsilonPMF[i][j][k][l], 0.1);	
					}
				}
			}
		}
		for(int l=0;l<computedLatLonMagEpsilonPMF[0][0][0].length;l++){
			System.out.println("Epsilon: "+(epsilonBinEdges[l]+epsilonBinEdges[l+1])/2);
			for(int k=0;k<computedLatLonMagEpsilonPMF[0][0].length;k++){
				System.out.println("Magnitude: "+(magBinEdges[k]+magBinEdges[k+1])/2);
				for(int j=0;j<computedLatLonMagEpsilonPMF[0].length;j++){
					System.out.print((lonBinEdges[j]+lonBinEdges[j+1])/2+" ");
				}
				System.out.println();
				for(int i=0;i<computedLatLonMagEpsilonPMF.length;i++){
					System.out.print((latBinEdges[i]+latBinEdges[i+1])/2+" ");
					for(int j=0;j<computedLatLonMagEpsilonPMF[i].length;j++){
						System.out.print(computedLatLonMagEpsilonPMF[i][j][k][l]+" ");
					}
					System.out.println();
				}
			}
		}
	}
	
	@Test
	public void getLatitudeLongitudeTectonicRegionTypePMFTest(){
		double[][][] computedLatLonTrtPMF = disCalc.getLatitudeLongitudeTectonicRegionTypePMF();
		double[][][] expectedLatLonTrtPMF = disCalcMC.getLatitudeLongitudeTectonicRegionTypePMF();
		
		for(int i=0;i<computedLatLonTrtPMF.length;i++){
			for(int j=0;j<computedLatLonTrtPMF[i].length;j++){
				for(int k=0;k<computedLatLonTrtPMF[i][j].length;k++){
					assertEquals(computedLatLonTrtPMF[i][j][k], expectedLatLonTrtPMF[i][j][k], 0.1);	
				}
			}
		}
	}
	
	@Test
	public void getLatitudeLongitudeMagnitudeEpsilonTectonicRegionTypePMFTest(){
		double[][][][][] computedDisMatrix = disCalc.getLatitudeLongitudeMagnitudeEpsilonTectonicRegionTypePMF();
		double[][][][][] expectedDisMatrix = disCalcMC.getLatitudeLongitudeMagnitudeEpsilonTectonicRegionTypePMF();
		
		for(int i=0;i<computedDisMatrix.length;i++){
			for(int j=0;j<computedDisMatrix[i].length;j++){
				for(int k=0;k<computedDisMatrix[i][j].length;k++){
					for(int l=0;l<computedDisMatrix[i][j][k].length;l++){
						for(int m=0;m<computedDisMatrix[i][j][k][l].length;m++){
							assertEquals(computedDisMatrix[i][j][k][l][m], expectedDisMatrix[i][j][k][l][m], 0.1);
						}
					}
				}
			}
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
