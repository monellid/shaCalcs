import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.GEM1ERF;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.util.TectonicRegionType;

/**
 * Disaggregation calculator based on Monte Carlo algorithm. Used for testing
 * classical disaggregation.
 * 
 * @author damianomonelli
 * 
 */
public class DisaggregationCalculatorMC {

	private double[] latBinEdges;
	private double[] lonBinEdges;
	private double[] magBinEdges;
	private double[] epsilonBinEdges;
	private double[] distanceBinEdges;
	private String[] tectonicRegionTypes;
	private double[][][][][] disaggregationMatrix;
	private Site site;

	// constructor
	public DisaggregationCalculatorMC(double[] latBinEdges,
			double[] lonBinEdges, double[] magBinEdges,
			double[] epsilonBinEdges, double[] distanceBinEdges) {
		this.latBinEdges = latBinEdges;
		this.lonBinEdges = lonBinEdges;
		this.magBinEdges = magBinEdges;
		this.epsilonBinEdges = epsilonBinEdges;
		this.distanceBinEdges = distanceBinEdges;
		this.tectonicRegionTypes = new String[] {
				TectonicRegionType.ACTIVE_SHALLOW.toString(),
				TectonicRegionType.STABLE_SHALLOW.toString(),
				TectonicRegionType.SUBDUCTION_INTERFACE.toString(),
				TectonicRegionType.SUBDUCTION_SLAB.toString(),
				TectonicRegionType.VOLCANIC.toString() };
		this.disaggregationMatrix = new double[latBinEdges.length - 1][lonBinEdges.length - 1][magBinEdges.length - 1][epsilonBinEdges.length - 1][tectonicRegionTypes.length];
	}

	public void disaggregate(
			double probExceed,
			Site site,
			GEM1ERF erf,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap,
			List<Double> imlVals, Random rn, int n) {

		this.site = site;

		// compute ground motion value corresponding to given probExceed
		double groundMotionValue;
		try {
			groundMotionValue = DisaggregationUtils.computeGroundMotionValue(
					probExceed, site, erf, imrMap, imlVals);
		} catch (Exception e) {
			throw new RuntimeException(e.toString());
		}

		// generate stochastic event set
		ArrayList<EqkRupture> ses = generateStochasticEventSets(erf, rn, n);

		// generate ground motion values for ruptures in ses, for the given site
		ArrayList<Double> gmfvs = computeGmfs(ses, imrMap, site, rn);

		int numRup = 0;
		for (int rupIndex = 0; rupIndex < ses.size();rupIndex++){
			
			EqkRupture rup = ses.get(rupIndex);

			// check if the rupture is inside the ranges considered
			// find closest point in the rupture area
			Location closestLoc = DisaggregationUtils.getClosestLocation(site,
					rup.getRuptureSurface());

			// get imr and set params
			ScalarIntensityMeasureRelationshipAPI imr = imrMap.get(rup
					.getTectRegType());
			imr.setSite(site);
			imr.setEqkRupture(rup);
			imr.setIntensityMeasureLevel(groundMotionValue);

			// rupture data
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

			// select
			// latitude-longitude-magnitude-tectonic_region_type-epsilon bin
			int ilat;
			for (ilat = 0; ilat < latBinEdges.length - 1; ilat++) {
				if (lat >= latBinEdges[ilat] && lat < latBinEdges[ilat + 1]) {
					break;
				}
			}
			int ilon;
			for (ilon = 0; ilon < lonBinEdges.length - 1; ilon++) {
				if (lon >= lonBinEdges[ilon] && lon < lonBinEdges[ilon + 1]) {
					break;
				}
			}
			int im;
			for (im = 0; im < magBinEdges.length - 1; im++) {
				if (magnitude >= magBinEdges[im]
						&& magnitude < magBinEdges[im + 1]) {
					break;
				}
			}
			int ie;
			for (ie = 0; ie < epsilonBinEdges.length - 1; ie++) {
				if (epsilon >= epsilonBinEdges[ie]
						&& epsilon < epsilonBinEdges[ie + 1]) {
					break;
				}
			}
			int itrt;
			for (itrt = 0; itrt < tectonicRegionTypes.length; itrt++) {
				if (tectonicRegionType
						.equalsIgnoreCase(tectonicRegionTypes[itrt])) {
					break;
				}
			}

			if (gmfvs.get(rupIndex) > groundMotionValue) {
				disaggregationMatrix[ilat][ilon][im][ie][itrt] = disaggregationMatrix[ilat][ilon][im][ie][itrt] + 1;
				numRup = numRup + 1;
			}

		}

		// normalize
		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
								disaggregationMatrix[i][j][k][l][m] = disaggregationMatrix[i][j][k][l][m]
										/numRup ;
						}
					}
				}
			}
		}

	}

	public double[] getMagnitudePMF() {

		double[] magPFM = new double[magBinEdges.length - 1];

		for (int i = 0; i < magBinEdges.length - 1; i++) {
			for (int j = 0; j < latBinEdges.length - 1; j++) {
				for (int k = 0; k < lonBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							magPFM[i] = magPFM[i]
									+ disaggregationMatrix[j][k][i][l][m];
						}
					}
				}
			}
		}

		return magPFM;
	}

	public double[] getDistancePMF() {
		double[] distPFM = new double[distanceBinEdges.length - 1];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							double meanLat = (latBinEdges[i] + latBinEdges[i + 1]) / 2;
							double meanLon = (lonBinEdges[j] + lonBinEdges[j + 1]) / 2;
							double distance = LocationUtils.horzDistance(site
									.getLocation(), new Location(meanLat,
									meanLon));
							if (distance < distanceBinEdges[0]
									|| distance > distanceBinEdges[distanceBinEdges.length - 1]) {
								continue;
							} else {
								int ii = 0;
								for (ii = 0; ii < distanceBinEdges.length - 1; ii++) {
									if (distance >= distanceBinEdges[ii]
											&& distance < distanceBinEdges[ii + 1])
										break;
								}
								distPFM[ii] = distPFM[ii]
										+ disaggregationMatrix[i][j][k][l][m];
							}
						}
					}
				}
			}
		}

		return distPFM;
	}

	public double[] getTectonicRegionTypePMF() {
		double[] trtPFM = new double[tectonicRegionTypes.length];

		for (int i = 0; i < tectonicRegionTypes.length; i++) {
			for (int j = 0; j < latBinEdges.length - 1; j++) {
				for (int k = 0; k < lonBinEdges.length - 1; k++) {
					for (int l = 0; l < magBinEdges.length - 1; l++) {
						for (int m = 0; m < epsilonBinEdges.length - 1; m++) {
							trtPFM[i] = trtPFM[i]
									+ disaggregationMatrix[j][k][l][m][i];
						}
					}
				}
			}
		}
		return trtPFM;
	}

	public double[][] getMagnitudeTectonicRegionTypePMF() {
		double[][] magTrtPMF = new double[magBinEdges.length - 1][tectonicRegionTypes.length];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							magTrtPMF[k][m] = magTrtPMF[k][m]
									+ disaggregationMatrix[i][j][k][l][m];
						}
					}
				}
			}
		}
		return magTrtPMF;
	}

	public double[][] getMagnitudeDistancePMF() {
		double[][] magDistPMF = new double[magBinEdges.length - 1][distanceBinEdges.length - 1];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							double meanLat = (latBinEdges[i] + latBinEdges[i + 1]) / 2;
							double meanLon = (lonBinEdges[j] + lonBinEdges[j + 1]) / 2;
							double distance = LocationUtils.horzDistance(site
									.getLocation(), new Location(meanLat,
									meanLon));
							if (distance < distanceBinEdges[0]
									|| distance > distanceBinEdges[distanceBinEdges.length - 1]) {
								continue;
							} else {
								int ii = 0;
								for (ii = 0; ii < distanceBinEdges.length - 1; ii++) {
									if (distance >= distanceBinEdges[ii]
											&& distance < distanceBinEdges[ii + 1])
										break;
								}
								magDistPMF[k][ii] = magDistPMF[k][ii]
										+ disaggregationMatrix[i][j][k][l][m];
							}
						}
					}
				}
			}
		}

		return magDistPMF;
	}

	public double[][][] getMagnitudeDistanceEpsilonPMF() {
		double[][][] magDistEpsilonPMF = new double[magBinEdges.length - 1][distanceBinEdges.length - 1][epsilonBinEdges.length - 1];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							double meanLat = (latBinEdges[i] + latBinEdges[i + 1]) / 2;
							double meanLon = (lonBinEdges[j] + lonBinEdges[j + 1]) / 2;
							double distance = LocationUtils.horzDistance(site
									.getLocation(), new Location(meanLat,
									meanLon));
							if (distance < distanceBinEdges[0]
									|| distance > distanceBinEdges[distanceBinEdges.length - 1]) {
								continue;
							} else {
								int ii = 0;
								for (ii = 0; ii < distanceBinEdges.length - 1; ii++) {
									if (distance >= distanceBinEdges[ii]
											&& distance < distanceBinEdges[ii + 1])
										break;
								}
								magDistEpsilonPMF[k][ii][l] = magDistEpsilonPMF[k][ii][l]
										+ disaggregationMatrix[i][j][k][l][m];
							}
						}
					}
				}
			}
		}
		return magDistEpsilonPMF;
	}

	public double[][] getLatitudeLongitudePMF() {
		double[][] latLonPMF = new double[latBinEdges.length - 1][lonBinEdges.length - 1];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							latLonPMF[i][j] = latLonPMF[i][j]
									+ disaggregationMatrix[i][j][k][l][m];
						}
					}
				}
			}
		}
		return latLonPMF;
	}

	public double[][][] getLatitudeLongitudeMagnitudePMF() {
		double[][][] latLonMagPMF = new double[latBinEdges.length - 1][lonBinEdges.length - 1][magBinEdges.length - 1];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							latLonMagPMF[i][j][k] = latLonMagPMF[i][j][k]
									+ disaggregationMatrix[i][j][k][l][m];
						}
					}
				}
			}
		}
		return latLonMagPMF;
	}

	public double[][][] getLatitudeLongitudeEpsilonPMF() {
		double[][][] latLonEpsilonPMF = new double[latBinEdges.length - 1][lonBinEdges.length - 1][epsilonBinEdges.length - 1];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							latLonEpsilonPMF[i][j][l] = latLonEpsilonPMF[i][j][l]
									+ disaggregationMatrix[i][j][k][l][m];
						}
					}
				}
			}
		}
		return latLonEpsilonPMF;
	}

	public double[][][][] getLatitudeLongitudeMagnitudeEpsilonPMF() {
		double[][][][] latLonMagEpsilonPMF = new double[latBinEdges.length - 1][lonBinEdges.length - 1][magBinEdges.length - 1][epsilonBinEdges.length - 1];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							latLonMagEpsilonPMF[i][j][k][l] = latLonMagEpsilonPMF[i][j][k][l]
									+ disaggregationMatrix[i][j][k][l][m];
						}
					}
				}
			}
		}
		return latLonMagEpsilonPMF;
	}

	public double[][][] getLatitudeLongitudeTectonicRegionTypePMF() {
		double[][][] latLonTrtPMF = new double[latBinEdges.length - 1][lonBinEdges.length - 1][tectonicRegionTypes.length];

		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							latLonTrtPMF[i][j][m] = latLonTrtPMF[i][j][m]
									+ disaggregationMatrix[i][j][k][l][m];
						}
					}
				}
			}
		}
		return latLonTrtPMF;
	}

	public double[][][][][] getLatitudeLongitudeMagnitudeEpsilonTectonicRegionTypePMF() {
		return disaggregationMatrix;
	}


	// generate stochastic event set
	private ArrayList<EqkRupture> generateStochasticEventSets(GEM1ERF erf,
			Random rn, int n) {
		ArrayList<EqkRupture> ses = new ArrayList<EqkRupture>();
		// generate n stochastic event sets
		for (int i = 0; i < n; i++) {
			ses.addAll(StochasticEventSetGenerator
					.getStochasticEventSetFromPoissonianERF(erf, rn));
		}
		return ses;
	}

	// generate ground motion values for a list of events on a site
	private ArrayList<Double> computeGmfs(
			ArrayList<EqkRupture> rupList,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap,
			Site site, Random rn) {

		ArrayList<Double> gmfs = new ArrayList<Double>();

		List<Site> sites = new ArrayList<Site>();
		sites.add(site);

		for (EqkRupture rup : rupList) {

			// compute ground motion field
			ScalarIntensityMeasureRelationshipAPI attenRel = imrMap.get(rup
					.getTectRegType());
			GroundMotionFieldCalculator gmfCalc = new GroundMotionFieldCalculator(
					attenRel, rup, sites);
			Map<Site, Double> gmf = gmfCalc
					.getUncorrelatedGroundMotionField(rn);
			gmfs.add(gmf.get(site));
		}
		return gmfs;
	}

}
