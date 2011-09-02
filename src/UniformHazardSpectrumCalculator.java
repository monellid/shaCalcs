import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFuncAPI;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.TectonicRegionType;

/**
 * Class implementing Uniform Hazard Spectrum (UHS) calculator. A UHS is defined
 * as spectral acceleration amplitudes, for a given number of periods,
 * corresponding to a given probability of exceedence.
 * 
 * @author damianomonelli
 * 
 */
public class UniformHazardSpectrumCalculator {

	private double[] periods;
	private EqkRupForecastAPI erf;
	private Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap;
	private List<Double> imlVals;

	/**
	 * Construct a UHS calculator for a specified list of periods.The erf,
	 * imrMap, and imlVals are utilized to computed the hazard curves from which
	 * the spectral values for the periods of interest are extracted. The ERF
	 * must be poissonian, if not an exception is thrown. Also checks that the
	 * maximum period requested for UHS is inside the range of periods supported
	 * by the considered GMPEs.
	 */
	public UniformHazardSpectrumCalculator(
			double[] periods,
			EqkRupForecastAPI erf,
			Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap,
			List<Double> imlVals) {
		this.periods = periods;
		this.erf = erf;
		this.imrMap = imrMap;
		this.imlVals = imlVals;

		// TODO: check also that the erf is poissonian
		validateInput(periods[periods.length - 1]);
	}

	/**
	 * Returns a UHS for a site and probability of exceedence.
	 */
	public double[] getUHS(double poE, Site site) {

		// UHS
		double[] uhs = new double[periods.length];

		// map containing hazard curves for each period
		Map<Double, DiscretizedFuncAPI> hazCurvesPerPeriod = new HashMap<Double, DiscretizedFuncAPI>();
		for (int i = 0; i < periods.length; i++) {
			DiscretizedFuncAPI hazCurve = new ArbitrarilyDiscretizedFunc();
			for (int j = 0; j < imlVals.size(); j++) {
				hazCurve.set(imlVals.get(j), 1.0);
			}
			hazCurvesPerPeriod.put(periods[i], hazCurve);
		}

		// loop over sources. For each source, loop over ruptures.
		for (int is = 0; is < erf.getNumSources(); is++) {

			ProbEqkSource src = erf.getSource(is);
			TectonicRegionType tectonicRegionType = src.getTectonicRegionType();

			// extract imr to be used, and set params
			ScalarIntensityMeasureRelationshipAPI imr = imrMap
					.get(tectonicRegionType);
			if (imr == null) {
				throw new RuntimeException("GMPE not defined for "
						+ tectonicRegionType.toString()
						+ ". Not able to continue calculation.");
			}
			imr.setSite(site);

			// extract period list
			SA_Param sa = (SA_Param) imr.getParameter(SA_Param.NAME);
			PeriodParam periodParam = sa.getPeriodParam();
			ArrayList<Double> periodList = periodParam.getAllowedDoubles();

			for (int ir = 0; ir < src.getNumRuptures(); ir++) {

				ProbEqkRupture rup = src.getRupture(ir);
				imr.setEqkRupture(rup);

				double rupProb = rup.getProbability();

				// compute probabilities of exceedence for the periods of
				// interest
				for (int ip = 0; ip < periods.length; ip++) {

					double period = periods[ip];

					// if period==0.0, this is PGA
					if (period == 0.0) {
						computeHazardCurveForPGA(hazCurvesPerPeriod, imr,
								rupProb, period);
					}
					// in case the period is equal to one of the supported
					// periods
					else if (periodParam.isAllowed(period)) {
						computeHazardCurveForSA(hazCurvesPerPeriod, imr,
								rupProb, period);

					} // if it is between PGA and the first SA
					else if (period > 0.0 && period < periodList.get(0)) {
						computeHazardCurveForPeriodBetweenPGAandSA(
								hazCurvesPerPeriod, imr, periodList, rupProb,
								period);

					} else {
						computeHazardCurveForInterpolatedPeriod(
								hazCurvesPerPeriod, imr, periodList, rupProb,
								period);

					}
				}

			}

		}

		return uhs;
	}

	private void computeHazardCurveForInterpolatedPeriod(
			Map<Double, DiscretizedFuncAPI> hazCurvesPerPeriod,
			ScalarIntensityMeasureRelationshipAPI imr,
			ArrayList<Double> periodList, double rupProb, double period) {
		
		// find the two closest periods and compute the
		// probability of exceedence by interpolating the
		// probabilities of exceedence of the closest periods
		int ii = getPeriodIndex(periodList, period);
		double period1 = periodList.get(ii);
		double period2 = periodList.get(ii + 1);
		imr.setIntensityMeasure(SA_Param.NAME);

		// loop over intensity measure levels
		for (int i = 0; i < imlVals.size(); i++) {
			imr.setIntensityMeasureLevel(imlVals.get(i));

			// compute probabilities of exceedence for the two
			// closest periods
			imr.getParameter(PeriodParam.NAME)
					.setValue(period1);
			double probExc1 = imr.getExceedProbability();
			imr.getParameter(PeriodParam.NAME)
					.setValue(period2);
			double probExc2 = imr.getExceedProbability();

			double probExc = getInterpolatedProb(period,
					period1, period2, probExc1, probExc2);

			computeHazardCurve(hazCurvesPerPeriod, rupProb, period, i, probExc);
		}
	}

	private void computeHazardCurve(
			Map<Double, DiscretizedFuncAPI> hazCurvesPerPeriod, double rupProb,
			double period, int i, double probExc) {
		double hazCurveValue = hazCurvesPerPeriod.get(
				period).getY(i);
		hazCurvesPerPeriod.get(period).set(
				i,
				hazCurveValue
						* Math.pow(1 - rupProb, probExc));
	}

	private int getPeriodIndex(ArrayList<Double> periodList, double period) {
		int ii;
		for (ii = 0; ii < periodList.size(); ii++) {
			if (period > periodList.get(ii)
					&& period < periodList.get(ii + 1)) {
				break;
			}
		}
		return ii;
	}

	private void computeHazardCurveForPeriodBetweenPGAandSA(
			Map<Double, DiscretizedFuncAPI> hazCurvesPerPeriod,
			ScalarIntensityMeasureRelationshipAPI imr,
			ArrayList<Double> periodList, double rupProb, double period) {

		double period1 = 0.0;
		double period2 = periodList.get(0);

		// loop over intensity measure levels
		for (int i = 0; i < imlVals.size(); i++) {
			imr.setIntensityMeasureLevel(imlVals.get(i));

			// compute probabilities of exceedence for PGA and first SA period
			imr.setIntensityMeasure(PGA_Param.NAME);
			double probExc1 = imr.getExceedProbability();
			imr.setIntensityMeasure(SA_Param.NAME);
			imr.getParameter(PeriodParam.NAME).setValue(period2);
			double probExc2 = imr.getExceedProbability();

			double probExc = getInterpolatedProb(period, period1, period2,
					probExc1, probExc2);

			computeHazardCurve(hazCurvesPerPeriod, rupProb, period, i, probExc);
		}
	}

	private void computeHazardCurveForSA(
			Map<Double, DiscretizedFuncAPI> hazCurvesPerPeriod,
			ScalarIntensityMeasureRelationshipAPI imr, double rupProb,
			double period) {

		imr.setIntensityMeasure(SA_Param.NAME);
		imr.getParameter(PeriodParam.NAME).setValue(period);

		// loop over intensity measure levels
		for (int i = 0; i < imlVals.size(); i++) {
			imr.setIntensityMeasureLevel(imlVals.get(i));
			double probExc = imr.getExceedProbability();
			computeHazardCurve(hazCurvesPerPeriod, rupProb, period, i, probExc);
		}
	}

	private void computeHazardCurveForPGA(
			Map<Double, DiscretizedFuncAPI> hazCurvesPerPeriod,
			ScalarIntensityMeasureRelationshipAPI imr, double rupProb,
			double period) {

		imr.setIntensityMeasure(PGA_Param.NAME);

		// loop over intensity measure levels
		for (int i = 0; i < imlVals.size(); i++) {
			imr.setIntensityMeasureLevel(imlVals.get(i));
			double probExc = imr.getExceedProbability();
			computeHazardCurve(hazCurvesPerPeriod, rupProb, period, i, probExc);
		}
	}

	// compute poE for period, given poEs for period1 and period2. Use linear
	// interpolation equation
	private double getInterpolatedProb(double period, double period1,
			double period2, double probExc1, double probExc2) {
		// interpolate
		// compute interpolation coefficients. Follow
		// equations
		// 3.3.1 and 3.3.2, pag. 107 in "Numerical Recipes
		// in
		// Fortran 77"
		double A = (period2 - period) / (period2 - period1);
		double B = 1 - A;
		double probExc = A * probExc1 + B * probExc2;
		return probExc;
	}

	// check that the specified GMPEs provide spectral accelerations
	// that cover the range of periods for UHS calculation
	private void validateInput(double maxPeriod) {

		for (ScalarIntensityMeasureRelationshipAPI imr : imrMap.values()) {

			// get list of supported periods
			SA_Param sa = (SA_Param) imr.getParameter(SA_Param.NAME);
			PeriodParam period = sa.getPeriodParam();
			ArrayList<Double> periodList = period.getAllowedDoubles();

			double supportedMaxPeriod = periodList.get(periodList.size() - 1);
			if (maxPeriod > supportedMaxPeriod) {
				throw new RuntimeException(
						"Maximum period requested for UHS is greater then maximum period supported by "
								+ imr.getName() + ". Not able to compute UHS.");
			}
		}

	}
}
