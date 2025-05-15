package lipid;

import adduct.Adduct;
import adduct.AdductList;
import lipid.Lipid;
import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IonizationMode ionizationMode;
    private String adduct; // !!TODO The adduct will be detected based on the groupedSignals
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;


    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IonizationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IonizationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        // !!TODO This set should be sorted according to help the program to deisotope the signals plus detect the adduct
        //In order for the TreeSet to work, the Peak class must implement Comparable, the compareTo method must be implemented to compare the mz values of the peaks
        //TreeSet uses the natural ordering of the elements — which is now defined in Peak.compareTo(...) as mz ascending.
        //Removes duplicates based on comparison (compareTo) not on equals()
        this.groupedSignals = new TreeSet<>(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IonizationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }


    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    // !CHECK Take into account that the score should be normalized between -1 and 1
    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    /**
     * @return The normalized score between 0 and 1 that consists on the final number divided into the times that the rule
     * has been applied.
     */
    public double getNormalizedScore() {
        return (double) this.score / this.totalScoresApplied;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }


    // !!TODO Detect the adduct with an algorithm or with drools, up to the user.
    public String detectAdductFromSignals(IonizationMode ionizationMode, double mzTolerance) {
        // Define la tolerancia en ppm para considerar dos masas como equivalentes


        if (groupedSignals == null || groupedSignals.size() < 2) {
            System.out.println("detectAdduct: Not enough signals (" +
                    (groupedSignals == null ? 0 : groupedSignals.size()) + ")");
            return "Unknown";
        }

        // Solo se ejecuta si la ionización es positiva
        if (ionizationMode == IonizationMode.POSITIVE) {

            // Recorremos todos los posibles aductos para la anotación principal
            for (Map.Entry<String, Double> entry : AdductList.MAPMZPOSITIVEADDUCTS.entrySet()) {

                String adductAnnotation = entry.getKey();
                double monoisotopicMassAdduct = Adduct.getMonoisotopicMassFromMZ(this.mz, adductAnnotation);

                System.out.println("\n--- Probando aducto para la anotación: " + adductAnnotation + " ---");
                System.out.println("Masa monoisotópica estimada de la anotación con " + adductAnnotation + ": " + monoisotopicMassAdduct);


                // Recorremos cada uno de los picos agrupados dentro de la anotación

                for (Peak peak : groupedSignals) {
                    double peakMz = peak.getMz();
                    System.out.println("\n   Analizando peak con m/z = " + peakMz);

                    if (Math.abs(peakMz - this.mz) <= mzTolerance) {
                        //It's the same peak as the objective, so I skip it.
                        continue;
                    }

                    // Probamos todos los posibles aductos para ese peak
                    for (Map.Entry<String, Double> entry2 : AdductList.MAPMZPOSITIVEADDUCTS.entrySet()) {
                        String adductPeak = entry2.getKey();
                        double mzForThatPeak = Adduct.getMZFromMonoisotopicMass(monoisotopicMassAdduct, adductPeak);
                        double diff = Math.abs(mzForThatPeak - peakMz);
                        System.out.println("secondAdduct=" + adductPeak + ", expectedMz=" + mzForThatPeak + ", observed=" + peakMz + ", diff=" + diff);
                        if (diff <= mzTolerance) {
                            System.out.println("DETECTED adduct: " + adductAnnotation + " (via " + adductPeak + ")");
                            return adductAnnotation;
                        }
                    }
                }


            }
        } else if (ionizationMode == IonizationMode.NEGATIVE) {
            // Recorremos todos los posibles aductos para la anotación principal
            for (Map.Entry<String, Double> entry : AdductList.MAPMZNEGATIVEADDUCTS.entrySet()) {

                String adductAnnotation = entry.getKey();
                double monoisotopicMassAdduct = Adduct.getMonoisotopicMassFromMZ(this.mz, adductAnnotation);

                System.out.println("\n--- Probando aducto para la anotación: " + adductAnnotation + " ---");
                System.out.println("Masa monoisotópica estimada de la anotación con " + adductAnnotation + ": " + monoisotopicMassAdduct);


                // Recorremos cada uno de los picos agrupados dentro de la anotación

                for (Peak peak : groupedSignals) {
                    double peakMz = peak.getMz();
                    System.out.println("\n   Analizando peak con m/z = " + peakMz);

                    if (Math.abs(peakMz - this.mz) <= mzTolerance) {
                        //It's the same peak as the objective, so I skip it.
                        continue;
                    }

                    // Probamos todos los posibles aductos para ese peak
                    for (Map.Entry<String, Double> entry2 : AdductList.MAPMZNEGATIVEADDUCTS.entrySet()) {
                        String adductPeak = entry2.getKey();
                        double mzForThatPeak = Adduct.getMZFromMonoisotopicMass(monoisotopicMassAdduct, adductPeak);
                        double diff = Math.abs(mzForThatPeak - peakMz);
                        System.out.println("secondAdduct=" + adductPeak + ", expectedMz=" + mzForThatPeak + ", observed=" + peakMz + ", diff=" + diff);
                        if (diff <= mzTolerance) {
                            System.out.println("DETECTED adduct: " + adductAnnotation + " (via " + adductPeak + ")");
                            return adductAnnotation;
                        }
                    }
                }
            }
        }
        return "Unkown";
    }

//    public String detectAdductFromGroupedSignals() {
//        double mzTolerance = 0.2;
//        if (groupedSignals == null || groupedSignals.size() < 2) {
//            System.out.println("detectAdduct: Not enough signals (" +
//                    (groupedSignals == null ? 0 : groupedSignals.size()) + ")");
//            return null;
//        }
//
//        Map<String, Double> adductMap = (getIonizationMode() == IonizationMode.POSITIVE)
//                ? AdductList.MAPMZPOSITIVEADDUCTS
//                : AdductList.MAPMZNEGATIVEADDUCTS;
//
//        double observedMz = this.getMz();
//        System.out.println("detectAdduct: observedMz = " + observedMz + ", mode = " + getIonizationMode());
//
//        // Para cada aducto candidato que podría explicar observedMz:
//        for (String candidateAdduct : adductMap.keySet()) {
//            System.out.println("  Probando candidateAdduct: " + candidateAdduct);
//            try {
//                // Masa monoisotópica según el aducto candidato para observedMz
//                double monoisotopicMass = Adduct.getMonoisotopicMassFromMZ(observedMz, candidateAdduct);
//                System.out.println("    monoisotopicMass = " + monoisotopicMass);
//
//                // Buscamos que alguno de los otros peaks se corresponda con esa misma masa
//                for (Peak otherPeak : groupedSignals) {
//                    System.out.println("    Comparando con otherPeak: " + otherPeak);
//                    if (Math.abs(otherPeak.getMz() - observedMz) <= mzTolerance) {
//                        // Es el mismo pico objetivo, lo saltamos
//                        System.out.println("      Skip: same as observedMz");
//                        continue;
//                    }
//                    // Probamos todos los aductos posibles para ese otherPeak
//                    for (String secondAdduct : adductMap.keySet()) {
//                        double expectedMz = Adduct.getMZFromMonoisotopicMass(monoisotopicMass, secondAdduct);
//                        double diff = Math.abs(expectedMz - otherPeak.getMz());
//                        System.out.println("      secondAdduct=" + secondAdduct + ", expectedMz=" + expectedMz + ", observed=" + otherPeak.getMz() + ", diff=" + diff);
//                        if (diff <= mzTolerance) {
//                            System.out.println("    DETECTED adduct: " + candidateAdduct + " (via " + secondAdduct + ")");
//                            return candidateAdduct;
//                        }
//                    }
//                }
//
//            } catch (IllegalArgumentException e) {
//                // Si el aducto no se puede parsear, simplemente lo ignoramos
//                System.out.println("    Ignorado candidateAdduct (parse error): " + candidateAdduct);
//            }
//        }
//
//        System.out.println("detectAdduct: Ningún aducto reconocido");
//        return null;
//    }
//
//


//    public static void main(String[] args) {
//        Peak mH = new Peak(700.500, 80000.0); // [M+H]+
//        Peak mNa = new Peak(722.482, 100000.0);  // [M+Na]+
//        Lipid lipid = new Lipid(1, "PC 34:1", "C42H82NO8P", "PC", 34, 1);
//
//        double annotationMZ = 722.482d;
//        double annotationIntensity = 100000.0;
//        double annotationRT = 6.5d;
//        Annotation annotation = new Annotation(lipid, annotationMZ, annotationIntensity, annotationRT, IonizationMode.POSITIVE, Set.of(mH, mNa));
//        System.out.println(annotation);
//
//        System.out.println("_");
//
//        Peak singlyCharged = new Peak(700.500, 100000.0);  // [M+H]+
//        Peak doublyCharged = new Peak(350.754, 85000.0);   // [M+2H]2+
//
//        Lipid lipid1 = new Lipid(3, "TG 54:3", "C57H104O6", "TG", 54, 3);
//        Annotation annotation1 = new Annotation(lipid1, singlyCharged.getMz(), singlyCharged.getIntensity(), 10d, IonizationMode.POSITIVE, Set.of(singlyCharged, doublyCharged));
//        System.out.println(annotation1);
//
//        System.out.println("_");
//
//        Peak mH2 = new Peak(700.500, 100000.0); // [M+H]+
//        Peak mk2 = new Peak(738.4564, 80000.0);  // [M+k]+
//        Lipid lipid2 = new Lipid(1, "PC 34:1", "C42H82NO8P", "PC", 34, 1);
//
//        double annotationRT2 = 6.5d;
//        Annotation annotation2 = new Annotation(lipid, mk2.getMz(), mk2.getIntensity(), annotationRT2, IonizationMode.POSITIVE, Set.of(mH, mk2));
//        System.out.println(annotation2);
//    }

}









