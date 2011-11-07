import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.param.DoubleParameter;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.griddedForecast.MagFreqDistsForFocalMechs;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.GEM1ERF;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.SourceData.GEMAreaSourceData;
import org.opensha.sha.earthquake.rupForecastImpl.GEM1.SourceData.GEMSourceData;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.attenRelImpl.BA_2008_AttenRel;
import org.opensha.sha.imr.attenRelImpl.CB_2008_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.ComponentParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

public class CalculatorsTestHelper {

	public static List<Double> getTestImlVals() {

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

	public static ScalarIntensityMeasureRelationshipAPI getTestIMRActiveShallow() {
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
	
	public static ScalarIntensityMeasureRelationshipAPI getTestIMRStableContinental() {
		ScalarIntensityMeasureRelationshipAPI imr = new CB_2008_AttenRel(null);
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

	public static double[] getTestIMRPeriods() {
		SA_Param sa = (SA_Param) CalculatorsTestHelper.getTestIMRActiveShallow()
				.getParameter(SA_Param.NAME);
		PeriodParam periodParam = sa.getPeriodParam();
		ArrayList<Double> periodValues = periodParam.getAllowedDoubles();
		double[] periods = new double[periodValues.size()];
		for (int i = 0; i < periods.length; i++) {
			periods[i] = periodValues.get(i);
		}
		return periods;
	}

	public static Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> getTestImrMap() {
		Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap = new HashMap<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI>();
		imrMap.put(TectonicRegionType.ACTIVE_SHALLOW, getTestIMRActiveShallow());
		imrMap.put(TectonicRegionType.STABLE_SHALLOW, getTestIMRStableContinental());
		return imrMap;
	}

	public static GEM1ERF getTestERF() {

		ArrayList<GEMSourceData> srcList = new ArrayList<GEMSourceData>();
		srcList.add(getTestSourceDataActiveShallow());
		srcList.add(getTestSourceDataStableCrust());
		double timeSpan = 50.0;
		//GEM1ERF erf = GEM1ERF.getGEM1ERF(srcList, timeSpan);
		//erf.getParameter(GEM1ERF.AREA_SRC_DISCR_PARAM_NAME).setValue(0.01);
		//erf.updateForecast();
		return GEM1ERF.getGEM1ERF(srcList, timeSpan);
	}

	public static GEMSourceData getTestSourceDataActiveShallow() {

		String id = "src1";
		String name = "testSource";
		TectonicRegionType tectReg = TectonicRegionType.ACTIVE_SHALLOW;
		LocationList border = new LocationList();
		border.add(new Location(-0.5, -0.5));
		border.add(new Location(0.5, -0.5));
		border.add(new Location(0.5, 0.5));
		border.add(new Location(-0.5, 0.5));
		Region reg = new Region(border, BorderType.MERCATOR_LINEAR);
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
		System.out.println("delta: "+magDist.getDelta());
		System.out.println("a value: "+magDist.get_aValue());
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
	
	public static GEMSourceData getTestSourceDataStableCrust() {

		String id = "src1";
		String name = "testSource";
		TectonicRegionType tectReg = TectonicRegionType.STABLE_SHALLOW;
		LocationList border = new LocationList();
		border.add(new Location(-0.5, 0.5));
		border.add(new Location(0.5, 0.5));
		border.add(new Location(0.5, 1.5));
		border.add(new Location(-0.5, 1.5));
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
	
	
	
	public static Site getTestSite(){
		Site site = new Site(new Location(0.0, 0.0));
		site.addParameter(new DoubleParameter(Vs30_Param.NAME, 760.0));
		site.addParameter(new DoubleParameter(DepthTo2pt5kmPerSecParam.NAME, 1.0));
		return site;
	}
	


}
