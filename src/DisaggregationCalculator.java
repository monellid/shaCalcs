import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.util.TectonicRegionType;

/**
 * Class implementing 'Disaggregation' calculator. Given an earthquake rupture
 * forecast {@link EqkRupForecastAPI} and a map of ground motion prediction
 * equations {@link ScalarIntensityMeasureRelationshipAPI}, it disaggregates the
 * ground motion value at a site corresponding to a given probability of
 * exceedance in a given time span. The calculator disaggregates in terms of
 * 
 * - latitude and longitude of the surface projection of the closest point of
 * the rupture area
 * 
 * - magnitude
 * 
 * - epsilon
 * 
 * - tectonic region type (active shallow crust, stable tectonic, subduction
 * interface, subduction intraslab, volcanic)
 * 
 * The class provides the method disaggregate() which is responsible for
 * calculating the 5 dimensional disaggregation matrix using a classical
 * approach, and disaggregateMonteCarlo(), that computes the same disaggregation
 * matrix but using a Monte Carlo approach. Additional methods are:
 * 
 * - getMagnitudePMF(): returns the 1D marginal probability mass function (PMF)
 * for magnitude
 * 
 * - getDistancePMF(): returns the 1D PMF for distance (where distance is the
 * closest distance to the surface projection of the rupture area)
 * 
 * -getTectonicRegionTypePMF(): returns the 1D PMF for tectonic region type
 * 
 * -getMagnitudeTectonicRegionTypePMF(): returns the 2D PMF for magnitude and
 * tectonic region type
 * 
 * - getMagnitudeDistancePMF(): returns the 2D PMF for magnitude and distance
 * 
 * - getMagnitudeDistanceEpsilonPMF(): returns the 3D PMF for magnitude,
 * distance, epsilon.
 * 
 * - getLatitudeLongitudePMF(): return the 2D PMF for latitude and longitude
 * 
 * - getLatitudeLongitudeMagnitudePMF(): return the 3D PMF for latitude and
 * longitude and magnitude
 * 
 * - getLatitudeLongitudeEpsilonPMF(): return the 3D PMF for latitude and
 * longitude and epsilon
 * 
 * - getLatitudeLongitudeMagnitudeEpsilonPMF(): return the 4D PMF for latitude,
 * longitude, magnitude and epsilon
 * 
 * - getLatitudeLongitudeTectonicRegionTypePMF(): return the 3D PMF for
 * latitude, longitude and tectonic region type
 * 
 * - getDisaggregationMatrix(): return the full 5D disaggregation matrix.
 * 
 * @author damianomonelli
 * 
 */
public class DisaggregationCalculator {

	private double[] latBinEdges;
	private double[] lonBinEdges;
	private double[] magBinEdges;
	private double[] epsilonBinEdges;
	private double[] distanceBinEdges;
	private String[] tectonicRegionTypes;
	private double[][][][][] disaggregationMatrix;
	private double minMag;

	/**
	 * The constructor accepts a list of latitudes, longitudes, magnitude, and
	 * epsilon values. The list of values are assumed to represent bin edges, to
	 * be used to compute the disaggregation matrix. Tectonic region types are
	 * initialized from the values in {@link TectonicRegionType}. The minimum
	 * magnitude is the one used in the ERF creation.
	 */
	public DisaggregationCalculator(double[] latBinEdges, double[] lonBinEdges,
			double[] magBinEdges, double[] epsilonBinEdges,
			double[] distanceBinEdges, double minMag) {
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
		this.minMag = minMag;
	}

	/**
	 * Disaggregates the intensity measure level corresponding to a given
	 * probability of exceedance in a given geographical location (site) in
	 * terms of latitude, longitude, magnitude, epsilon, and tectonic region
	 * type. It uses a classical approach: that is contribution are computed
	 * directly from probabilities. The earthquake rupture forecast is required
	 * to be Poissonian. If not an exception is thrown. An exception is also
	 * thrown if a GMPEs is set to have zero standard deviation (in this case
	 * epsilon cannot be computed). Returns thre ground motion value
	 * corresponding to the given probability of exceedence.
	 */
	public double disaggregate(
			double probExceed,
			Site site,
			EqkRupForecastAPI erf,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap,
			List<Double> imlVals) {

		CalculatorsUtils.ensurePoissonian(erf);

		CalculatorsUtils.ensureNonZeroStd(imrMap);

		// compute ground motion value corresponding to given probExceed
		double groundMotionValue = CalculatorsUtils.computeGroundMotionValue(
				probExceed, site, erf, imrMap, imlVals);

		double totalAnnualRate = 0.0;

		// loop over sources. For each source compute annual rate of ruptures
		// with magnitude greater then Mmin (lambda(M>Mmin). For each source,
		// loop over ruptures. For each rupture Rup_i compute probability of
		// exceeding the given IML (that is P(IML > x, Rup_i)). The annual rate
		// of events producing an IML > x, with given magnitude, latitude,
		// longitude, epsilon and tectonic region type, is then computed as:
		// lambda(IML>x, M=m, LAT=lat, LON=lon, EPSILON=epsilon,
		// TECTONIC_REGIONTYPE=tectonicRegionType) = lambda(M>Mmin) * P(IML > x,
		// Rup_i) * P(Rup_i)
		for (int i = 0; i < erf.getNumSources(); i++) {

			ProbEqkSource src = erf.getSource(i);

			// compute total rate above Mmin
			double totProb = src.computeTotalProbAbove(minMag);
			double totRate = -Math.log(1 - totProb);

			// select IMR based on source tectonic region type
			ScalarIntensityMeasureRelationshipAPI imr = imrMap.get(src
					.getTectonicRegionType());
			imr.setSite(site);
			imr.setIntensityMeasureLevel(groundMotionValue);

			// loop over ruptures
			for (int j = 0; j < src.getNumRuptures(); j++) {

				ProbEqkRupture rup = src.getRupture(j);
				imr.setEqkRupture(rup);

				// find closest point in the rupture area
				Location closestLoc = CalculatorsUtils.getClosestLocation(site,
						rup.getRuptureSurface());

				// get rupture parameters
				double lat = closestLoc.getLatitude();
				double lon = closestLoc.getLongitude();
				double magnitude = rup.getMag();
				double epsilon = imr.getEpsilon();
				String tectonicRegionType = src.getTectonicRegionType()
						.toString();

				// if one of the parameters (latitude, longitude, magnitude,
				// epsilon) is outside of
				// the considered range do not include in the conditional
				// probability calculation, that is skip the rupture
				if (!isInsideRange(lat, lon, magnitude, epsilon)) {
					continue;
				}

				// select
				// latitude-longitude-magnitude-tectonic_region_type-epsilon bin
				int[] binIndex = getBinIndex(lat, lon, magnitude, epsilon,
						tectonicRegionType);

				// compute annual rate of events exceeding the given iml
				double annualRate = totRate * imr.getExceedProbability()
						* rup.getProbability();

				// update disaggregation matrix
				disaggregationMatrix[binIndex[0]][binIndex[1]][binIndex[2]][binIndex[3]][binIndex[4]] = disaggregationMatrix[binIndex[0]][binIndex[1]][binIndex[2]][binIndex[3]][binIndex[4]]
						+ annualRate;

				// sum rates to get total annual rate
				totalAnnualRate = totalAnnualRate + annualRate;
			}

		}

		// normalize by total annual rate
		normalizaDisaggregationMatrix(totalAnnualRate);

		return groundMotionValue;
	}

	/**
	 * Compute disaggregation matrix using Monte Carlo approach. That is
	 * generates multiple (n) stochastic event sets, and for each rupture
	 * compute a stochastic realization of the ground motion field in the site
	 * of interest.
	 */
	public double disaggregateMonteCarlo(
			double probExceed,
			Site site,
			EqkRupForecastAPI erf,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap,
			List<Double> imlVals, Random rn, int n) {

		// compute ground motion value corresponding to given probExceed
		double groundMotionValue = CalculatorsUtils.computeGroundMotionValue(
				probExceed, site, erf, imrMap, imlVals);

		// generate stochastic event set
		ArrayList<EqkRupture> ses = generateStochasticEventSets(erf, rn, n);

		// generate ground motion values for ruptures in ses, for the given site
		ArrayList<Double> gmfvs = computeGmfs(ses, imrMap, site, rn);

		int numRup = 0;
		for (int rupIndex = 0; rupIndex < ses.size(); rupIndex++) {

			EqkRupture rup = ses.get(rupIndex);

			// check if the rupture is inside the ranges considered
			// find closest point in the rupture area
			Location closestLoc = CalculatorsUtils.getClosestLocation(site,
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
			if (!isInsideRange(lat, lon, magnitude, epsilon)) {
				continue;
			}

			// select
			// latitude-longitude-magnitude-tectonic_region_type-epsilon bin
			int[] binIndex = getBinIndex(lat, lon, magnitude, epsilon,
					tectonicRegionType);

			if (gmfvs.get(rupIndex) > groundMotionValue) {
				disaggregationMatrix[binIndex[0]][binIndex[1]][binIndex[2]][binIndex[3]][binIndex[4]] = disaggregationMatrix[binIndex[0]][binIndex[1]][binIndex[2]][binIndex[3]][binIndex[4]] + 1;
				numRup = numRup + 1;
			}

		}

		// normalize by total number of ruptures exceeding the expected ground
		// motion value
		normalizaDisaggregationMatrix(numRup);

		return groundMotionValue;
	}

	private void normalizaDisaggregationMatrix(double normalizationFactor) {
		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							disaggregationMatrix[i][j][k][l][m] = disaggregationMatrix[i][j][k][l][m]
									/ normalizationFactor;
						}
					}
				}
			}
		}
	}

	private boolean isInsideRange(double lat, double lon, double magnitude,
			double epsilon) {
		if (lat < latBinEdges[0] || lat >= latBinEdges[latBinEdges.length - 1]) {
			return false;
		}
		if (lon < lonBinEdges[0] || lon >= lonBinEdges[lonBinEdges.length - 1]) {
			return false;
		}
		if (magnitude < magBinEdges[0]
				|| magnitude >= magBinEdges[magBinEdges.length - 1]) {
			return false;
		}
		if (epsilon < epsilonBinEdges[0]
				|| epsilon >= epsilonBinEdges[epsilonBinEdges.length - 1]) {
			return false;
		}
		return true;
	}

	private int[] getBinIndex(double lat, double lon, double magnitude,
			double epsilon, String tectonicRegionType) {
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
			if (magnitude >= magBinEdges[im] && magnitude < magBinEdges[im + 1]) {
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
			if (tectonicRegionType.equalsIgnoreCase(tectonicRegionTypes[itrt])) {
				break;
			}
		}

		return new int[] { ilat, ilon, im, ie, itrt };
	}

	/**
	 * Extract magnitude PMF from disaggregation matrix.
	 */
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

	/**
	 * Extract distance PMF from disaggregation matrix. Distance is calculated
	 * with respect to a site.
	 */
	public double[] getDistancePMF(Site site) {
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

	/**
	 * Extract tectonic region type PMF from disaggregation matrix. The tectonic
	 * region types are those define in {@link TectonicRegionType}.
	 */
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

	/**
	 * Extract magnitude-tectonic region type PMF from disaggregation matrix.
	 */
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

	/**
	 * Extract magnitude-distance PMF from disaggregation matrix. Distance is
	 * calculated with respect to a site.
	 */
	public double[][] getMagnitudeDistancePMF(Site site) {
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

	/**
	 * Extract magnitude-distance-epsilon PMF from disaggregation matrix.
	 * Distance is calculated with respect to a site.
	 */
	public double[][][] getMagnitudeDistanceEpsilonPMF(Site site) {
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

	/**
	 * Extract latitude-longitude PMF from disaggregation matrix.
	 */
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

	/**
	 * Extract latitude-longitude-magnitude PMF from disaggregation matrix.
	 */
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

	/**
	 * Extract latitude-longitude-epsilon PMF from disaggregation matrix.
	 */
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

	/**
	 * Extract latitude-longitude-magnitude-epsilon PMF from disaggregation
	 * matrix.
	 */
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

	/**
	 * Extract latitude-longitude-tectonic region type PMF from disaggregation
	 * matrix.
	 */
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

	/**
	 * Returns full disaggregation matrix
	 * (latitude-longitude-magnitude-epsilon-tectonic region type).
	 */
	public double[][][][][] getDisaggregationMatrix() {
		return disaggregationMatrix;
	}

	/**
	 * Generate multiple (n) stochastic event sets.
	 */
	private ArrayList<EqkRupture> generateStochasticEventSets(
			EqkRupForecastAPI erf, Random rn, int n) {
		ArrayList<EqkRupture> ses = new ArrayList<EqkRupture>();
		// generate n stochastic event sets
		for (int i = 0; i < n; i++) {
			ses.addAll(StochasticEventSetGenerator
					.getStochasticEventSetFromPoissonianERF(erf, rn));
		}
		return ses;
	}

	/**
	 * Generate stochastic ground motion values for a list of events on a site.
	 */
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
