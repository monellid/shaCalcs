import java.util.ArrayList;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.util.TectonicRegionType;

public class DisaggregationCalculator {

	private double[] magBinEdges;
	private double[] distBinEdges;
	private double[] epsilonBinEdges;

	/**
	 * The constructor accepts a list of magnitude, distance and epsilon values.
	 * Values are assumed to represent bin centers for computing the conditional
	 * probability distribution of M, R, and epsilon. M, R, and epsilon are
	 * assumed to be equally spaced.
	 * TODO: change to accept bin edges (more flexible)
	 */
	public DisaggregationCalculator(double[] mag, double[] dist,
			double[] epsilon) {
		// TODO: check that values are equally spaced
		this.magBinEdges = new double[mag.length + 1];
		this.distBinEdges = new double[dist.length + 1];
		this.epsilonBinEdges = new double[epsilon.length + 1];
		setBinEdges(mag, this.magBinEdges);
		setBinEdges(dist, this.distBinEdges);
		setBinEdges(epsilon, this.epsilonBinEdges);
	}

	// compute bin edges from bin center values
	private void setBinEdges(double[] val, double[] binEdges) {
		double deltaMag = val[1] - val[0];
		for (int i = 0; i < binEdges.length - 1; i++) {
			binEdges[i] = val[i] - deltaMag / 2;
			binEdges[i + 1] = val[i] + deltaMag / 2;
		}
	}

	/**
	 * Disaggregates an intensity measure level (iml, input as ln(iml) when
	 * using a GMPE or simply iml when using an IPE) in a given geographical
	 * location (site) in terms of magnitude, distance, epsilon. The earthquake
	 * rupture forecast is assumed to be Poissonian. Distance is the closest
	 * distance to the rupture.
	 */
	public double[][][] getMagDistEpsilonDisaggregation(
			double iml,
			Site site,
			EqkRupForecastAPI erf,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap) {

		// ensurePoissonian(erf);

		// TODO: make sure that non-zero standard deviation is set in imr

		DistanceRupParameter distRup = new DistanceRupParameter();

		double[][][] conditionalPMF = new double[magBinEdges.length - 1][distBinEdges.length - 1][epsilonBinEdges.length - 1];
		for (int i = 0; i < magBinEdges.length - 1; i++) {
			for (int j = 0; j < distBinEdges.length - 1; j++) {
				for (int k = 0; k < epsilonBinEdges.length - 1; k++) {
					conditionalPMF[i][j][k] = 0.0;
				}
			}
		}

		double totalAnnualRate = 0.0;

		// loop over sources. For each source, loop over ruptures. For each
		// rupture, compute the mean annual rate of exceedance (the specified
		// iml), and accumulate it in the corresponding
		// magnitude-distance-epsilon bin.
		// The mean annual rate of exceedance is computed as:
		// annualRate = -Math.log(1.0 - probExceedance)
		// where probExceedance is the probability of exceedance (of the
		// specified
		// iml), given the rupture, and site characteristic. The equation comes
		// from the equation for poissonian probability: P(x>=iml) = 1 - exp(-nu
		// * T), solved for nu and assuming T = 1.
		// Finally normalize with the total mean annual rate of exceedance
		// as given by all the ruptures. The total mean annual rate is obtained
		// by summing the mean annual rates associated to each rupture.
		for (int i = 0; i < erf.getNumSources(); i++) {

			ProbEqkSource src = erf.getSource(i);

			// select IMR based on source tectonic region type
			ScalarIntensityMeasureRelationshipAPI imr = imrMap.get(src
					.getTectonicRegionType());
			imr.setSite(site);
			imr.setIntensityMeasureLevel(iml);

			// loop over ruptures
			for (int j = 0; j < src.getNumRuptures(); j++) {

				ProbEqkRupture rup = src.getRupture(j);
				imr.setEqkRupture(rup);

				// if magnitude, distance and epsilon are outside of the
				// considered range do not include in the conditional
				// probability calculation
				double magnitude = rup.getMag();
				double distance = (Double) distRup.getValue(rup, site);
				double epsilon = imr.getEpsilon();
				if (magnitude < magBinEdges[0]
						|| magnitude >= magBinEdges[magBinEdges.length - 1]) {
					continue;
				}
				if (distance < distBinEdges[0]
						|| distance >= distBinEdges[distBinEdges.length - 1]) {
					continue;
				}
				if (epsilon < epsilonBinEdges[0]
						|| epsilon >= epsilonBinEdges[epsilonBinEdges.length - 1]) {
					continue;
				}

				// select magnitude-distance-epsilon bin
				int im;
				for (im = 0; im < magBinEdges.length - 1; im++) {
					if (magnitude >= magBinEdges[im]
							&& magnitude < magBinEdges[im + 1]) {
						break;
					}
				}
				int id;
				for (id = 0; id < distBinEdges.length - 1; id++) {
					if (distance >= distBinEdges[id]
							&& distance < distBinEdges[id + 1]) {
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

				double probExceedance = imr.getExceedProbability();
				double annualRate = -Math.log(1.0 - probExceedance);

				conditionalPMF[im][id][ie] = conditionalPMF[im][id][ie]
						+ annualRate;

				totalAnnualRate = totalAnnualRate + annualRate;
			}
		}

		// normalize by total annual rate
		for (int i = 0; i < magBinEdges.length - 1; i++) {
			for (int j = 0; j < distBinEdges.length - 1; j++) {
				for (int k = 0; k < epsilonBinEdges.length - 1; k++) {
					conditionalPMF[i][j][k] = conditionalPMF[i][j][k]
							/ totalAnnualRate;
				}
			}
		}

		return conditionalPMF;
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
