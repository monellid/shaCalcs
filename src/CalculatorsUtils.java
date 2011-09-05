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
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.EvenlyGriddedSurfaceAPI;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.util.TectonicRegionType;

public class CalculatorsUtils {

	/**
	 * Ensure ERF contains only Poissonian sources
	 */
	public static void ensurePoissonian(EqkRupForecastAPI erf) {
		for (ProbEqkSource src : (ArrayList<ProbEqkSource>) erf.getSourceList())
			if (src.isSourcePoissonian() == false)
				throw new IllegalArgumentException("Sources must be Poissonian");
	}

	/**
	 * Ensure non-zero standard deviation in
	 * {@link ScalarIntensityMeasureRelationshipAPI} map
	 */
	public static void ensureNonZeroStd(
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap) {
		for (ScalarIntensityMeasureRelationshipAPI imr : imrMap.values()) {
			String stdDevType = (String) imr.getParameter(StdDevTypeParam.NAME)
					.getValue();
			if (stdDevType.equalsIgnoreCase(StdDevTypeParam.STD_DEV_TYPE_NONE)) {
				throw new RuntimeException(
						"Attenuation relationship must have a non-zero standard deviation");
			}
		}
	}

	/**
	 * compute ground motion value corresponding to probability of exceedance on
	 * a site, given an earthquake rupture forecast, a GMPE map, and a list of
	 * IML values.
	 */
	public static double computeGroundMotionValue(
			double probExceed,
			Site site,
			EqkRupForecastAPI erf,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap,
			List<Double> imlVals) {
		DiscretizedFuncAPI hazardCurve = new ArbitrarilyDiscretizedFunc();
		for (double val : imlVals)
			hazardCurve.set(val, 1.0);
		HazardCurveCalculator hcc;
		try {
			hcc = new HazardCurveCalculator();
			hcc.getHazardCurve(hazardCurve, site, imrMap, erf);
		} catch (Exception e) {
			throw new RuntimeException(e.toString());
		}
		return getGroundMotionValueForPoE(probExceed, hazardCurve);
	}

	/**
	 * Extract ground motion value corresponding to a probability of Exceedence
	 * from hazard curve.
	 */
	public static double getGroundMotionValueForPoE(double probExceed,
			DiscretizedFuncAPI hazardCurve) {
		double groundMotionValue;
		if (probExceed > hazardCurve.getY(0)) {
			groundMotionValue = hazardCurve.getX(0);
		} else if (probExceed < hazardCurve.getY(hazardCurve.getNum() - 1)) {
			groundMotionValue = hazardCurve.getX(hazardCurve.getNum() - 1);
		} else {
			groundMotionValue = hazardCurve.getFirstInterpolatedX(probExceed);
		}
		return groundMotionValue;
	}

	/** Compute rupture's closest location to site.
	 */
	public static Location getClosestLocation(Site site,
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

}
