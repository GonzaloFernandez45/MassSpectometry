package lipid;

import org.junit.Before;
import org.junit.Test;

import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

public class AdductDetectionTest {
    // !!TODO For the adduct detection both regular algorithms or drools can be used as far the tests are passed.


    @Before
    public void setup() {
        // !! TODO Empty by now,you can create common objects for all tests.
    }

    @Test
    public void shouldDetectAdductBasedOnMzDifference() {

        // Given two peaks with ~21.98 Da difference (e.g., [M+H]+ and [M+Na]+)
        Peak mH = new Peak(700.500, 100000.0); // [M+H]+
        Peak mNa = new Peak(722.482, 80000.0);  // [M+Na]+
        Lipid lipid = new Lipid(1, "PC 34:1", "C42H82NO8P", "PC", 34, 1);

        double annotationMZ = 700.49999d;
        double annotationIntensity = 80000.0;
        double annotationRT = 6.5d;
        Annotation annotation = new Annotation(lipid, annotationMZ, annotationIntensity, annotationRT, IonizationMode.POSITIVE, Set.of(mH, mNa));


        // Then we should call the algorithmic/knowledge system rules fired to detect the adduct and Set it!
        //
         String adduct = annotation.detectAdductFromSignals(IonizationMode.POSITIVE, 10d);
        annotation.setAdduct(adduct);

        assertNotNull("[M+H]+ should be detected", annotation.getAdduct());
        assertEquals( "Adduct inferred from lowest mz in group","[M+H]+", annotation.getAdduct());
    }


    @Test
    public void shouldDetectLossOfWaterAdduct() {
        Peak mh = new Peak(700.500, 90000.0);        // [M+H]+
        Peak mhH2O = new Peak(682.4894, 70000.0);     // [M+H–H₂O]+, ~18.0106 Da less

        Lipid lipid = new Lipid(1, "PE 36:2", "C41H78NO8P", "PE", 36, 2);
        Annotation annotation = new Annotation(lipid, mh.getMz(), mh.getIntensity(), 7.5d, IonizationMode.POSITIVE, Set.of(mh, mhH2O));

        String adduct = annotation.detectAdductFromSignals(IonizationMode.POSITIVE, 10d);
        annotation.setAdduct(adduct);


        assertNotNull("[M+H]+ should be detected", annotation.getAdduct());

        assertEquals( "Adduct inferred from lowest mz in group","[M+H]+", annotation.getAdduct());
    }

    @Test
    public void shouldDetectDoublyChargedAdduct() {
        // Assume real M = (700.500 - 1.0073) = 699.4927
        // So [M+2H]2+ = (M + 2.0146) / 2 = 350.7536
        Peak singlyCharged = new Peak(701.500, 100000.0);  // [M+H]+
        Peak doublyCharged = new Peak(350.754, 85000.0);   // [M+2H]2+

        Lipid lipid = new Lipid(3, "TG 54:3", "C57H104O6", "TG", 54, 3);
        Annotation annotation = new Annotation(lipid, singlyCharged.getMz(), singlyCharged.getIntensity(), 10d, IonizationMode.POSITIVE, Set.of(singlyCharged, doublyCharged));

        String adduct = annotation.detectAdductFromSignals(IonizationMode.POSITIVE,10d);
        annotation.setAdduct(adduct);

        assertNotNull("[M+H]+ should be detected", annotation.getAdduct());

        assertEquals( "Adduct inferred from lowest mz in group","[M+H]+", annotation.getAdduct());
    }

    @Test
    public void shouldDetectAdductBasedOnMzk() {

        // Given two peaks with ~21.98 Da difference (e.g., [M+H]+ and [M+Na]+)
        Peak mH = new Peak(700.500, 100000.0); // [M+H]+
        Peak mk = new Peak(738.4564, 80000.0);  // [M+k]+
        Lipid lipid = new Lipid(1, "PC 34:1", "C42H82NO8P", "PC", 34, 1);

        double annotationRT = 6.5d;
        Annotation annotation = new Annotation(lipid, mk.getMz(), mk.getIntensity(), annotationRT, IonizationMode.POSITIVE, Set.of(mH, mk));

        String adduct = annotation.detectAdductFromSignals(IonizationMode.POSITIVE, 2d);
        annotation.setAdduct(adduct);

        assertNotNull("[M+K]+ should be detected", annotation.getAdduct());
        assertEquals( "Adduct inferred from lowest mz in group","[M+K]+", annotation.getAdduct());
        System.out.println(annotation.getAdduct());

    }

    @Test
    public void shouldDetectLossOfWaterAdduct1() {
        Peak mh = new Peak(700.500, 90000.0);        // [M+H]+
        Peak mhH2O = new Peak(682.4894, 70000.0);     // [M+H–H₂O]+, ~18.0106 Da less

        Lipid lipid = new Lipid(1, "PE 36:2", "C41H78NO8P", "PE", 36, 2);
        Annotation annotation = new Annotation(lipid, mhH2O.getMz(), mhH2O.getIntensity(), 7.5d, IonizationMode.POSITIVE, Set.of(mh, mhH2O));

        String adduct = annotation.detectAdductFromSignals(IonizationMode.POSITIVE, 0.2d);
        annotation.setAdduct(adduct);

        assertNotNull("[M+H-H2O]+ should be detected", annotation.getAdduct());

        assertEquals( "Adduct inferred from lowest mz in group","[M+H-H2O]+", annotation.getAdduct());
    }

    @Test
    public void shouldDetectDoublyChargedAdduct2() {
        // Assume real M = (700.500 - 1.0073) = 699.4927
        // So [M+2H]2+ = (M + 2.0146) / 2 = 350.7536
        Peak singlyCharged = new Peak(701.500, 100000.0);  // [M+H]+
        Peak doublyCharged = new Peak(350.754, 85000.0);   // [M+2H]2+

        Lipid lipid = new Lipid(3, "TG 54:3", "C57H104O6", "TG", 54, 3);
        Annotation annotation = new Annotation(lipid, doublyCharged.getMz(), doublyCharged.getIntensity(), 10d, IonizationMode.POSITIVE, Set.of(singlyCharged, doublyCharged));


        String adduct = annotation.detectAdductFromSignals(IonizationMode.POSITIVE,10d);
        annotation.setAdduct(adduct);

        assertNotNull("[M+2H]2+ should be detected", annotation.getAdduct());

        assertEquals( "Adduct inferred from lowest mz in group","[M+2H]2+", annotation.getAdduct());
    }


    //NEGATIVE IONIZATION MODE TESTS

    @Test
    public void shouldDetectNegativeAdductBasedOnMzDifference() {
        // [M–H]– y [M+Cl]– tienen una diferencia de ~36.96 Da
        Peak mH = new Peak(700.500, 100000.0);      // [M–H]–
        Peak mCl = new Peak(737.460, 80000.0);      // [M+Cl]–

        Lipid lipid = new Lipid(4, "PI 38:4", "C47H83O13P", "PI", 38, 4);
        Annotation annotation = new Annotation(lipid, mH.getMz(), mH.getIntensity(), 8.0, IonizationMode.NEGATIVE, Set.of(mH, mCl));

        String adduct = annotation.detectAdductFromSignals(IonizationMode.NEGATIVE, 10d);
        annotation.setAdduct(adduct);

        assertNotNull("[M-H]− should be detected", annotation.getAdduct());
        assertEquals("[M-H]−", annotation.getAdduct());
    }

    @Test
    public void shouldDetectNegativeLossOfWaterAdduct() {
        // [M–H]– y [M–H–H2O]– diferencia de ~18.0106 Da
        Peak mh = new Peak(700.500, 90000.0);         // [M–H]–
        Peak mhH2O = new Peak(682.4894, 70000.0);      // [M–H–H₂O]–

        Lipid lipid = new Lipid(5, "PS 36:2", "C42H76NO10P", "PS", 36, 2);
        Annotation annotation = new Annotation(lipid, mhH2O.getMz(), mhH2O.getIntensity(), 7.5, IonizationMode.NEGATIVE, Set.of(mh, mhH2O));

        String adduct = annotation.detectAdductFromSignals(IonizationMode.NEGATIVE, 5d);
        annotation.setAdduct(adduct);

        assertNotNull("[M-H-H2O]− should be detected", annotation.getAdduct());
        assertEquals("[M-H-H2O]−", annotation.getAdduct());
    }

    @Test
    public void shouldDetectNegativeAdductWithFormicAcid() {
        // [M–H]– y [M+HCOOH-H]– diferencia de ~44.9982 Da
        Peak mH = new Peak(600.300, 95000.0);              // [M–H]–
        Peak mFormic = new Peak(645.2982, 85000.0);        // [M+HCOOH-H]– = 600.300 + 44.9982

        Lipid lipid = new Lipid(6, "PE 36:4", "C41H72NO8P", "PE", 36, 4);
        Annotation annotation = new Annotation(lipid, mFormic.getMz(), mFormic.getIntensity(), 6.8, IonizationMode.NEGATIVE, Set.of(mH, mFormic));

        String adduct = annotation.detectAdductFromSignals(IonizationMode.NEGATIVE, 5d);
        annotation.setAdduct(adduct);

        assertNotNull("[M+HCOOH-H]− should be detected", annotation.getAdduct());
        assertEquals("[M+HCOOH-H]−", annotation.getAdduct());
    }


}
