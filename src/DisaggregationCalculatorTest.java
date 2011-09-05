import static org.junit.Assert.*;

import java.util.List;
import java.util.Map;
import java.util.Random;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.data.Site;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.GEM1ERF;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
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
	private DisaggregationCalculator disCalcMC;

	@Before
	public void setUp() {

		latBinEdges = new double[] { -0.6, -0.3, -0.1, 0.1, 0.3, 0.6 };
		lonBinEdges = new double[] { -0.6, -0.3, -0.1, 0.1, 0.3, 0.6 };
		magBinEdges = new double[] { 5.0, 6.0, 7.0, 8.0, 9.0 };
		epsilonBinEdges = new double[] { -3.5, -2.5, -1.5, -0.5, +0.5, +1.5,
				+2.5, +3.5 };
		distanceBinEdges = new double[] { 0, 20, 40, 60 };

		site = CalculatorsTestHelper.getTestSite();
		erf = CalculatorsTestHelper.getTestERF();
		imrMap = CalculatorsTestHelper.getTestImrMap();
		imlVals = CalculatorsTestHelper.getTestImlVals();

		double minMag = (Double) (erf.getParameter(GEM1ERF.MIN_MAG_NAME)
				.getValue());

		// compute disaggregation matrix following classical approach
		disCalc = new DisaggregationCalculator(latBinEdges, lonBinEdges,
				magBinEdges, epsilonBinEdges, distanceBinEdges, minMag);
		disCalc.disaggregate(probExceed, site, erf, imrMap, imlVals);

		// compute disaggregation calculator with MC approach
		disCalcMC = new DisaggregationCalculator(latBinEdges, lonBinEdges,
				magBinEdges, epsilonBinEdges, distanceBinEdges, minMag);
		disCalcMC.disaggregateMonteCarlo(probExceed, site, erf, imrMap,
				imlVals, rn, n);
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

		double[] computedDistancePMF = disCalc.getDistancePMF(site);

		double[] expectedDistancePMF = disCalcMC.getDistancePMF(site);

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
			System.out.println("tectonic region type: " + (i + 1)
					+ ", computed: " + computedTectonicRegionTypePMF[i]
					+ ", expected: " + expectedTectonicRegionTypePMF[i]);
			assertEquals(computedTectonicRegionTypePMF[i],
					expectedTectonicRegionTypePMF[i], 0.1);
		}
	}

	@Test
	public void getMagnitudeTectonicRegionTypePMFTest() {
		double[][] computedMagTrtPMF = disCalc
				.getMagnitudeTectonicRegionTypePMF();

		double[][] expectedMagTrtPMF = disCalcMC
				.getMagnitudeTectonicRegionTypePMF();

		for (int i = 0; i < computedMagTrtPMF.length; i++) {
			for (int j = 0; j < computedMagTrtPMF[i].length; j++) {
				assertEquals(computedMagTrtPMF[i][j], expectedMagTrtPMF[i][j],
						0.1);
			}
		}

	}

	@Test
	public void getMagnitudeDistancePMFTest() {
		double[][] computedMagDistPMF = disCalc.getMagnitudeDistancePMF(site);

		double[][] expectedMagDistPMF = disCalcMC.getMagnitudeDistancePMF(site);

		for (int i = 0; i < computedMagDistPMF.length; i++) {
			for (int j = 0; j < computedMagDistPMF[i].length; j++) {
				assertEquals(computedMagDistPMF[i][j],
						expectedMagDistPMF[i][j], 0.1);
			}
		}

		for (int j = 0; j < computedMagDistPMF[0].length; j++) {
			System.out.print((distanceBinEdges[j] + distanceBinEdges[j + 1])
					/ 2 + " ");
		}
		System.out.println();
		for (int i = 0; i < computedMagDistPMF.length; i++) {
			System.out.print((magBinEdges[i] + magBinEdges[i + 1]) / 2 + " ");
			for (int j = 0; j < computedMagDistPMF[i].length; j++) {
				System.out.print(computedMagDistPMF[i][j] + " ");
			}
			System.out.println();
		}
	}

	@Test
	public void getMagnitudeDistanceEpsilonPMFTest() {
		double[][][] computedMagDistEpsilonPMF = disCalc
				.getMagnitudeDistanceEpsilonPMF(site);

		double[][][] expectedMagDistEpsilonPMF = disCalcMC
				.getMagnitudeDistanceEpsilonPMF(site);

		for (int i = 0; i < computedMagDistEpsilonPMF.length; i++) {
			for (int j = 0; j < computedMagDistEpsilonPMF[i].length; j++) {
				for (int k = 0; k < computedMagDistEpsilonPMF[i][j].length; k++) {
					assertEquals(computedMagDistEpsilonPMF[i][j][k],
							expectedMagDistEpsilonPMF[i][j][k], 0.1);
				}
			}
		}
		for (int k = 0; k < computedMagDistEpsilonPMF[0][0].length; k++) {
			System.out.println("Epsilon: "
					+ (epsilonBinEdges[k] + epsilonBinEdges[k + 1]) / 2);
			for (int j = 0; j < computedMagDistEpsilonPMF[0].length; j++) {
				System.out
						.print((distanceBinEdges[j] + distanceBinEdges[j + 1])
								/ 2 + " ");
			}
			System.out.println();
			for (int i = 0; i < computedMagDistEpsilonPMF.length; i++) {
				System.out.print((magBinEdges[i] + magBinEdges[i + 1]) / 2
						+ " ");
				for (int j = 0; j < computedMagDistEpsilonPMF[i].length; j++) {
					System.out.print(computedMagDistEpsilonPMF[i][j][k] + " ");
				}
				System.out.println();
			}
		}
	}

	@Test
	public void getLatituteLongitudePMFTest() {
		double[][] computedLatLonPMF = disCalc.getLatitudeLongitudePMF();
		double[][] expectedLatLonPMF = disCalcMC.getLatitudeLongitudePMF();

		for (int i = 0; i < computedLatLonPMF.length; i++) {
			for (int j = 0; j < computedLatLonPMF[i].length; j++) {
				assertEquals(computedLatLonPMF[i][j], expectedLatLonPMF[i][j],
						0.1);
			}
		}

		for (int j = 0; j < computedLatLonPMF[0].length; j++) {
			System.out.print((latBinEdges[j] + latBinEdges[j + 1]) / 2 + " ");
		}
		System.out.println();
		for (int i = 0; i < computedLatLonPMF.length; i++) {
			System.out.print((lonBinEdges[i] + lonBinEdges[i + 1]) / 2 + " ");
			for (int j = 0; j < computedLatLonPMF[i].length; j++) {
				System.out.print(computedLatLonPMF[i][j] + " ");
			}
			System.out.println();
		}
	}

	@Test
	public void getLatitudeLongitudeMagnitudePMFTest() {
		double[][][] computedLatLonMagPMF = disCalc
				.getLatitudeLongitudeMagnitudePMF();
		double[][][] expectedLatLonMagPMF = disCalcMC
				.getLatitudeLongitudeMagnitudePMF();

		for (int i = 0; i < computedLatLonMagPMF.length; i++) {
			for (int j = 0; j < computedLatLonMagPMF[i].length; j++) {
				for (int k = 0; k < computedLatLonMagPMF[i][j].length; k++) {
					assertEquals(computedLatLonMagPMF[i][j][k],
							expectedLatLonMagPMF[i][j][k], 0.1);
				}
			}
		}
		for (int k = 0; k < computedLatLonMagPMF[0][0].length; k++) {
			System.out.println("Magnitude: "
					+ (magBinEdges[k] + magBinEdges[k + 1]) / 2);
			for (int j = 0; j < computedLatLonMagPMF[0].length; j++) {
				System.out.print((lonBinEdges[j] + lonBinEdges[j + 1]) / 2
						+ " ");
			}
			System.out.println();
			for (int i = 0; i < computedLatLonMagPMF.length; i++) {
				System.out.print((latBinEdges[i] + latBinEdges[i + 1]) / 2
						+ " ");
				for (int j = 0; j < computedLatLonMagPMF[i].length; j++) {
					System.out.print(computedLatLonMagPMF[i][j][k] + " ");
				}
				System.out.println();
			}
		}
	}

	@Test
	public void getLatitudeLongitudeEpsilonPMFTest() {
		double[][][] computedLatLonEpsilonPMF = disCalc
				.getLatitudeLongitudeEpsilonPMF();
		double[][][] expectedLatLonEpsilonPMF = disCalcMC
				.getLatitudeLongitudeEpsilonPMF();

		for (int i = 0; i < computedLatLonEpsilonPMF.length; i++) {
			for (int j = 0; j < computedLatLonEpsilonPMF[i].length; j++) {
				for (int k = 0; k < computedLatLonEpsilonPMF[i][j].length; k++) {
					assertEquals(computedLatLonEpsilonPMF[i][j][k],
							expectedLatLonEpsilonPMF[i][j][k], 0.1);
				}
			}
		}
		for (int k = 0; k < computedLatLonEpsilonPMF[0][0].length; k++) {
			System.out.println("Epsilon: "
					+ (epsilonBinEdges[k] + epsilonBinEdges[k + 1]) / 2);
			for (int j = 0; j < computedLatLonEpsilonPMF[0].length; j++) {
				System.out.print((lonBinEdges[j] + lonBinEdges[j + 1]) / 2
						+ " ");
			}
			System.out.println();
			for (int i = 0; i < computedLatLonEpsilonPMF.length; i++) {
				System.out.print((latBinEdges[i] + latBinEdges[i + 1]) / 2
						+ " ");
				for (int j = 0; j < computedLatLonEpsilonPMF[i].length; j++) {
					System.out.print(computedLatLonEpsilonPMF[i][j][k] + " ");
				}
				System.out.println();
			}
		}
	}

	@Test
	public void getLatitudeLongitudeMagnitudeEpsilonPMFTest() {
		double[][][][] computedLatLonMagEpsilonPMF = disCalc
				.getLatitudeLongitudeMagnitudeEpsilonPMF();
		double[][][][] expectedLatLonMagEpsilonPMF = disCalcMC
				.getLatitudeLongitudeMagnitudeEpsilonPMF();

		for (int i = 0; i < computedLatLonMagEpsilonPMF.length; i++) {
			for (int j = 0; j < computedLatLonMagEpsilonPMF[i].length; j++) {
				for (int k = 0; k < computedLatLonMagEpsilonPMF[i][j].length; k++) {
					for (int l = 0; l < computedLatLonMagEpsilonPMF[i][j][k].length; l++) {
						assertEquals(computedLatLonMagEpsilonPMF[i][j][k][l],
								expectedLatLonMagEpsilonPMF[i][j][k][l], 0.1);
					}
				}
			}
		}
		for (int l = 0; l < computedLatLonMagEpsilonPMF[0][0][0].length; l++) {
			System.out.println("Epsilon: "
					+ (epsilonBinEdges[l] + epsilonBinEdges[l + 1]) / 2);
			for (int k = 0; k < computedLatLonMagEpsilonPMF[0][0].length; k++) {
				System.out.println("Magnitude: "
						+ (magBinEdges[k] + magBinEdges[k + 1]) / 2);
				for (int j = 0; j < computedLatLonMagEpsilonPMF[0].length; j++) {
					System.out.print((lonBinEdges[j] + lonBinEdges[j + 1]) / 2
							+ " ");
				}
				System.out.println();
				for (int i = 0; i < computedLatLonMagEpsilonPMF.length; i++) {
					System.out.print((latBinEdges[i] + latBinEdges[i + 1]) / 2
							+ " ");
					for (int j = 0; j < computedLatLonMagEpsilonPMF[i].length; j++) {
						System.out
								.print(computedLatLonMagEpsilonPMF[i][j][k][l]
										+ " ");
					}
					System.out.println();
				}
			}
		}
	}

	@Test
	public void getLatitudeLongitudeTectonicRegionTypePMFTest() {
		double[][][] computedLatLonTrtPMF = disCalc
				.getLatitudeLongitudeTectonicRegionTypePMF();
		double[][][] expectedLatLonTrtPMF = disCalcMC
				.getLatitudeLongitudeTectonicRegionTypePMF();

		for (int i = 0; i < computedLatLonTrtPMF.length; i++) {
			for (int j = 0; j < computedLatLonTrtPMF[i].length; j++) {
				for (int k = 0; k < computedLatLonTrtPMF[i][j].length; k++) {
					assertEquals(computedLatLonTrtPMF[i][j][k],
							expectedLatLonTrtPMF[i][j][k], 0.1);
				}
			}
		}
	}

	@Test
	public void getLatitudeLongitudeMagnitudeEpsilonTectonicRegionTypePMFTest() {
		double[][][][][] computedDisMatrix = disCalc
				.getDisaggregationMatrix();
		double[][][][][] expectedDisMatrix = disCalcMC
				.getDisaggregationMatrix();

		for (int i = 0; i < computedDisMatrix.length; i++) {
			for (int j = 0; j < computedDisMatrix[i].length; j++) {
				for (int k = 0; k < computedDisMatrix[i][j].length; k++) {
					for (int l = 0; l < computedDisMatrix[i][j][k].length; l++) {
						for (int m = 0; m < computedDisMatrix[i][j][k][l].length; m++) {
							assertEquals(computedDisMatrix[i][j][k][l][m],
									expectedDisMatrix[i][j][k][l][m], 0.1);
						}
					}
				}
			}
		}
	}
}
