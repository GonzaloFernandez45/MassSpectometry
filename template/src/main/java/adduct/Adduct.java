package adduct;

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

        double massToSearch;
        int multimer = extractMultimer(adduct); // regex "(\\d+)M"
        int charge = extractCharge(adduct);     // regex "(\\d+)([+-]$)"
        Double adductMass = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        if(adductMass == null) {
            adductMass = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
        }
        if(adductMass == null) {
            adductMass = 0.0;
        }
        massToSearch = (mz + adductMass) * charge / multimer;
        return massToSearch;
    }

    public static void main(String[] args) {
        // Test the methods
        Double mz = 700.500;
        String adduct = "[M+H]+";
        Double monoisotopicMass = getMonoisotopicMassFromMZ(mz, adduct);
        System.out.println("Monoisotopic mass: " + monoisotopicMass);

        Double calculatedMz = getMZFromMonoisotopicMass(monoisotopicMass, adduct);
        System.out.println("Calculated mz: " + calculatedMz);

        int ppmIncrement = calculatePPMIncrement(mz, monoisotopicMass);
        System.out.println("PPM Increment: " + ppmIncrement);

        double deltaPPM = calculateDeltaPPM(mz, 10);
        System.out.println("Delta PPM: " + deltaPPM);
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
            throw new IllegalArgumentException("Adduct not found: " + adduct);
        }
        if (charge == 1) {
            //if Adduct is single charge the formula is M = m/z +- adductMass. Charge is 1 so it does not affect
            // therefore mz = M + adductMass
            mz = monoisotopicMass + adductMass;
        } else if (multimer > 1) {
            //if adduct is a dimer or multimer the formula is M =  (mz +- adductMass) / numberOfMultimer
            //therefore mz = (M * numberOfMultimer) + adductMass
            mz = (monoisotopicMass*multimer) + adductMass ;
        } else {
            //if Adduct is double or triple charged the formula is M = ( mz +- adductMass ) * charge
            //therefore mz = (M / charge) + adductMass
            mz = (monoisotopicMass/charge) + adductMass;
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
    private static int extractMultimer(String adduct) {
        Pattern pMultimer = Pattern.compile("(\\d+)M");
        Matcher mMultimer = pMultimer.matcher(adduct);
        if (mMultimer.find()) {
            return Integer.parseInt(mMultimer.group(1));
        }else {
            return 1;
        }
    }
    private static int extractCharge(String adduct) {
        Pattern pCharge = Pattern.compile("(\\d+)([+-]$)");
        Matcher mCharge = pCharge.matcher(adduct);
        if (mCharge.find()) {
            return Integer.parseInt(mCharge.group(1));
        } else {
            return 1;
        }
    }




}
