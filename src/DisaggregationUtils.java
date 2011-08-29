import java.rmi.RemoteException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFuncAPI;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.faultSurface.EvenlyGriddedSurfaceAPI;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.util.TectonicRegionType;

public class DisaggregationUtils {

	// compute ground motion value corresponding to probability of exceedance
	public static double computeGroundMotionValue(
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

	// compute rupture closest location to site
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
