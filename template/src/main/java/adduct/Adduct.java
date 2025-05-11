package adduct;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Adduct {

    /*
    DATOS IMPROTANTES

        if Adduct is single charge the formula is M = m/z +- adductMass. Charge is 1 so it does not affect

        if Adduct is double or triple charged the formula is M = ( mz +- adductMass ) * charge

        if adduct is a dimer or multimer the formula is M =  (mz +- adductMass) / numberOfMultimer

        return monoisotopicMass;

         */

    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */

    // Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).

    public static Double getMonoisotopicMassFromMZ(Double mz, String adduct) {


        Double mzoffset = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        if(mzoffset == null) {
            mzoffset = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
        }
        if(mzoffset == null) {
            mzoffset = 0.0;
        }

        int multimer = extractMultimer(adduct); // regex "(\\d+)M"
        int charge = extractCharge(adduct);     // regex "(\\d+)([+-]$)"
         return (mz + mzoffset) * charge / multimer;
    }



    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoisotopicMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return mz
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct) {
        Double mz;

        int multimer = extractMultimer(adduct);
        int charge = extractCharge(adduct);

        Double adductMass = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        if (adductMass == null) {
            adductMass = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
        }

        if (adductMass == null) {
            adductMass = 0.0;
        }
        if (charge == 1) {
            //if Adduct is single charge the formula is M = m/z +- adductMass. Charge is 1 so it does not affect
            // therefore mz = M + adductMass
            mz = monoisotopicMass - adductMass;
        } else if (multimer > 1) {
            //if adduct is a dimer or multimer the formula is M =  (mz +- adductMass) / numberOfMultimer
            //therefore mz = (M * numberOfMultimer) + adductMass
            mz = (monoisotopicMass*multimer) - adductMass ;
        } else {
            //if Adduct is double or triple charged the formula is M = ( mz +- adductMass ) * charge
            //therefore mz = (M / charge) + adductMass
            mz = (monoisotopicMass/charge) - adductMass;
        }
        return mz;
    }





    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass) * 1000000
                / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param ppm ppm of tolerance
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        deltaPPM =  Math.round(Math.abs((experimentalMass * ppm) / 1000000));
        return deltaPPM;

    }

    // metodos extra que hacen falta:
    public static int extractMultimer(String adduct) {
        Pattern pMultimer = Pattern.compile("(\\d+)M");
        Matcher mMultimer = pMultimer.matcher(adduct);
        if (mMultimer.find()) {
            return Integer.parseInt(mMultimer.group(1));
        }else {
            return 1;
        }
    }
    public static int extractCharge(String adduct) {
        Pattern pCharge = Pattern.compile("(\\d+)([+-]$)");
        Matcher mCharge = pCharge.matcher(adduct);
        if (mCharge.find()) {
            return Integer.parseInt(mCharge.group(1));
        } else {
            return 1;
        }
    }

    public static void main(String[] args) {
        List<Double> monoisotopicmasses = new java.util.ArrayList<>();
        monoisotopicmasses.add(699.49);
        monoisotopicmasses.add(1396.98);
        monoisotopicmasses.add(677.5107);
        monoisotopicmasses.add(661.5368);
        monoisotopicmasses.add(682.4894);
        monoisotopicmasses.add(717.503);
        monoisotopicmasses.add(1362.9167);
        monoisotopicmasses.add(349.746);
        monoisotopicmasses.add(338.755);

        List<Double> ppms = new java.util.ArrayList<>();
        List<Double> ppms2 = new java.util.ArrayList<>();
        List<Double> ppms3 = new java.util.ArrayList<>();

        double mzPeak = 700.500;
        double mzPeak2 = mzPeak*2;
        double mzPeak3 = mzPeak/2;
        double monoisotopicMassPeak = getMonoisotopicMassFromMZ(mzPeak, "");
        double monoisotopicMassPeak2 = getMonoisotopicMassFromMZ(mzPeak2, "");
        double monoisotopicMassPeak3 = getMonoisotopicMassFromMZ(mzPeak3, "");
        System.out.println("Monoisotopic Mass Peak: " + monoisotopicMassPeak);
        System.out.println("Monoisotopic Mass Peak2: " + monoisotopicMassPeak2);
        System.out.println("Monoisotopic Mass Peak3: " + monoisotopicMassPeak3);


        for (Double monoisotopicMass : monoisotopicmasses) {
            double ppm = calculatePPMIncrement(monoisotopicMass, monoisotopicMassPeak);
            double ppm2 = calculatePPMIncrement(monoisotopicMass, monoisotopicMassPeak2);
            double ppm3 = calculatePPMIncrement(monoisotopicMass, monoisotopicMassPeak3);

            ppms.add(ppm);
            ppms2.add(ppm2);
            ppms3.add(ppm3);
        }
        System.out.println("PPM: " + ppms);
        System.out.println("PPM2: " + ppms2);
        System.out.println("PPM3: " + ppms3);
        System.out.println("\n");


        List<Double> _ppms = new java.util.ArrayList<>();
        List<Double> _ppms2 = new java.util.ArrayList<>();
        List<Double> _ppms3 = new java.util.ArrayList<>();

        List<Double> monoisotopicmasses2 = new java.util.ArrayList<>();
        monoisotopicmasses2.add(349.74);
        monoisotopicmasses2.add(697.47);
        monoisotopicmasses2.add(327.76);
        monoisotopicmasses2.add(311.79);
        monoisotopicmasses2.add(332.72);
        monoisotopicmasses2.add(367.757);
        monoisotopicmasses2.add(663.42);
        monoisotopicmasses2.add(174.873);
        monoisotopicmasses2.add(163.882);

        double _mzPeak = 350.754;
        double _mzPeak2 = _mzPeak*2;
        double _mzPeak3 = _mzPeak/2;
        double _monoisotopicMassPeak = getMonoisotopicMassFromMZ(_mzPeak, "");
        double _monoisotopicMassPeak2 = getMonoisotopicMassFromMZ(_mzPeak2, "");
        double _monoisotopicMassPeak3 = getMonoisotopicMassFromMZ(_mzPeak3, "");
        System.out.println("Monoisotopic Mass Peak: " + _monoisotopicMassPeak);
        System.out.println("Monoisotopic Mass Peak2: " + _monoisotopicMassPeak2);
        System.out.println("Monoisotopic Mass Peak3: " + _monoisotopicMassPeak3);


        for (Double monoisotopicMass : monoisotopicmasses2) {
            double _ppm = calculatePPMIncrement(monoisotopicMass, _monoisotopicMassPeak);
            double _ppm2 = calculatePPMIncrement(monoisotopicMass, _monoisotopicMassPeak2);
            double _ppm3 = calculatePPMIncrement(monoisotopicMass, _monoisotopicMassPeak3);

            _ppms.add(_ppm);
            _ppms2.add(_ppm2);
            _ppms3.add(_ppm3);
        }
        System.out.println("PPM: " + _ppms);
        System.out.println("PPM2: " + _ppms2);
        System.out.println("PPM3: " + _ppms3);

    }


}
