import java.rmi.RemoteException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFuncAPI;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.GEM1ERF;
import org.opensha.sha.faultSurface.EvenlyGriddedSurfaceAPI;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
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
 * calculating the 5 dimensional disaggregation matrix. Additional methods are:
 * 
 * - getMagnitudePMF(): returns the 1D marginal probability mass function (PMF)
 * for magnitude
 * 
 * - getDistancePMF(): returns the 1D marginal probability mass function (PMF)
 * for distance (where distance is the closest distance to the surface
 * projection of the rupture area)
 * 
 * -getTectonicRegionTypePMF(): returns the 1D marginal probability mass
 * function (PMF) for tectonic region type
 * 
 * - getMagnitudeDistancePMF(): returns the 2d marginal probability mass
 * function (PMF) for magnitude and distance (again closest distance to the
 * surface projection of the rupture area)
 * 
 * - getMagnitudeDistanceEpsilonPMF(): returns the 3d marginal probability mass
 * function (PMF) for magnitude, distance, epsilon.
 * 
 * - getLatitudeLongitudePMF(): return the 2d marginal probability mass function
 * (PMF) for latitude and longitude
 * 
 * - getMagnitudeTectonicRegionTypePMF(): return the 2d marginal probability
 * mass function (PMF) for magnitude and tectonic region type
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
	private Site site;

	/**
	 * The constructor accepts a list of latitudes, longitudes, magnitude, and
	 * epsilon values. The list of values are assumed to represent bin edges, to
	 * be used to compute the disaggregation matrix. Tectonic region types are
	 * initialized from the values in {@link TectonicRegionType}.
	 */
	public DisaggregationCalculator(double[] latBinEdges, double[] lonBinEdges,
			double[] magBinEdges, double[] epsilonBinEdges,
			double[] distanceBinEdges) {
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

	/**
	 * Disaggregates the intensity measure level corresponding to a given
	 * probability of exceedance in a given geographical location (site) in
	 * terms of latitude, longitude, magnitude, epsilon, and tectonic region
	 * type. The earthquake rupture forecast is assumed to be Poissonian.
	 * 
	 * @throws RemoteException
	 */
	public void disaggregate(
			double probExceed,
			Site site,
			GEM1ERF erf,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap,
			List<Double> imlVals) throws RemoteException {

		this.site = site;

		// ensure that the ERF contains only Poissonian sources
		ensurePoissonian(erf);

		// TODO: make sure that non-zero standard deviation is set in imr. If it
		// is set to zero, then the imr.getEpsilon() method returns infinity.

		// compute ground motion value corresponding to given probExceed
		double groundMotionValue = computeGroundMotionValue(probExceed, site,
				erf, imrMap, imlVals);

		// initialize disaggregation matrix
		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							disaggregationMatrix[i][j][k][l][m] = 0.0;
						}
					}
				}
			}
		}

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
			double totProb = src.computeTotalProbAbove((Double) (erf
					.getParameter(GEM1ERF.MIN_MAG_NAME).getValue()));
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
				Location closestLoc = getClosestLocation(site,
						rup.getRuptureSurface());

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

				double annualRate = totRate * imr.getExceedProbability()
						* rup.getProbability();

				disaggregationMatrix[ilat][ilon][im][ie][itrt] = disaggregationMatrix[ilat][ilon][im][ie][itrt]
						+ annualRate;

				totalAnnualRate = totalAnnualRate + annualRate;
			}

		}

		// normalize by total annual rate
		for (int i = 0; i < latBinEdges.length - 1; i++) {
			for (int j = 0; j < lonBinEdges.length - 1; j++) {
				for (int k = 0; k < magBinEdges.length - 1; k++) {
					for (int l = 0; l < epsilonBinEdges.length - 1; l++) {
						for (int m = 0; m < tectonicRegionTypes.length; m++) {
							disaggregationMatrix[i][j][k][l][m] = disaggregationMatrix[i][j][k][l][m]
									/ totalAnnualRate;
						}
					}
				}
			}
		}
	}

	public double[] getMagnitudePMF() {

		double[] magPFM = new double[magBinEdges.length - 1];
		for (int i = 0; i < magPFM.length; i++) {
			magPFM[i] = 0.0;
		}

		for (int i = 0; i < magPFM.length; i++) {
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
		for (int i = 0; i < distPFM.length; i++) {
			distPFM[i] = 0.0;
		}

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
		return trtPFM;
	}

	public Location getClosestLocation(Site site,
			EvenlyGriddedSurfaceAPI rupSurf) {
		Location closestLoc = null;
		double closestDistance = Double.MAX_VALUE;
		for (Location loc : rupSurf.getLocationList()) {
			double distance = Math
					.sqrt(Math.pow(
							LocationUtils.horzDistance(site.getLocation(), loc),
							2)
							+ Math.pow(
									LocationUtils.vertDistance(
											site.getLocation(), loc), 2));
			if (distance < closestDistance) {
				closestDistance = distance;
				closestLoc = loc;
			}
		}
		return closestLoc;
	}

	public double computeGroundMotionValue(
			double probExceed,
			Site site,
			EqkRupForecastAPI erf,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap,
			List<Double> imlVals) throws RemoteException {
		double groundMotionValue;
		DiscretizedFuncAPI hazardCurve = new ArbitrarilyDiscretizedFunc();
		for (double val : imlVals)
			hazardCurve.set(val, 1.0);
		HazardCurveCalculator hcc = new HazardCurveCalculator();
		hcc.getHazardCurve(hazardCurve, site, imrMap, erf);
		// this assumes that the hazard curves values (imlVals) are in ascending
		// order
		// TODO: check for this
		if (probExceed > hazardCurve.getY(0)) {
			groundMotionValue = hazardCurve.getX(0);
		} else if (probExceed < hazardCurve.getY(hazardCurve.getNum() - 1)) {
			groundMotionValue = hazardCurve.getX(hazardCurve.getNum() - 1);
		} else {
			groundMotionValue = hazardCurve.getFirstInterpolatedX(probExceed);
		}
		return groundMotionValue;
	}

	/**
	 * Check if the ERF contains only Poissonian sources
	 */
	private Boolean ensurePoissonian(EqkRupForecastAPI erf) {
		for (ProbEqkSource src : (ArrayList<ProbEqkSource>) erf.getSourceList())
			if (src.isSourcePoissonian() == false)
				throw new IllegalArgumentException("Sources must be Poissonian");
		return true;
	}
}
