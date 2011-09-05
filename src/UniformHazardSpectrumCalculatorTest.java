import java.util.List;
import java.util.Map;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.sha.earthquake.EqkRupForecastAPI;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.util.TectonicRegionType;

public class UniformHazardSpectrumCalculatorTest {

	private double[] periods;
	private EqkRupForecastAPI erf;
	private Map<TectonicRegionType, ScalarIntensityMeasureRelationshipAPI> imrMap;
	private List<Double> imlVals;
	private double maxDistance;
	private UniformHazardSpectrumCalculator uhsCalc;

	@Before
	public void setUp() {
		periods = CalculatorsTestHelper.getTestIMRPeriods();
		erf = CalculatorsTestHelper.getTestERF();
		imrMap = CalculatorsTestHelper.getTestImrMap();
		imlVals = CalculatorsTestHelper.getTestImlVals();
		maxDistance = 200.0;
		uhsCalc = new UniformHazardSpectrumCalculator(
				periods, erf, imrMap, imlVals, maxDistance);
	}
	
	@After
	public void tearDown(){
		periods = null;
		erf = null;
		imrMap = null;
		imlVals = null;
		uhsCalc = null;
	}
	
	@Test
	public void UHSTest(){
		double poE = 0.1;
		double[] uhs = uhsCalc.getUHS(poE, CalculatorsTestHelper.getTestSite());
		for(int i=0;i<periods.length;i++){
			System.out.println(periods[i]+" "+Math.exp(uhs[i]));
		}
	}

}
