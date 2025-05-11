package lipid;

import adduct.Adduct;
import adduct.AdductList;
import org.drools.base.rule.QueryArgument;

import java.sql.SQLOutput;
import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IoniationMode ionizationMode;
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
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode) {
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
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode, Set<Peak> groupedSignals) {
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

    public IoniationMode getIonizationMode() {
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
    public String detectAdductFromSignals(IoniationMode ionizationMode) {
        double tolerance = 2000; // ppm tolerance for adduct detection

        if (ionizationMode == IoniationMode.POSITIVE) {
            for (Peak peak : groupedSignals) {
                System.out.println("peak: " + peak.getMz());
                for (Map.Entry<String, Double> entry : AdductList.MAPMZPOSITIVEADDUCTS.entrySet()) {
                    double monoisotopicMassAnnotation = Adduct.getMonoisotopicMassFromMZ(mz, entry.getKey());


                    System.out.println("Adduct: " + entry.getKey());
                    System.out.println("mz: "+ mz);
                    System.out.println("Monoisotopic Mass Annotation: " + monoisotopicMassAnnotation);

                    double monoisotopicMassPeak = Adduct.getMonoisotopicMassFromMZ(peak.getMz(), entry.getKey());
                     System.out.println("Monoisotopic Mass Peak: " + monoisotopicMassPeak);
                    int ppmIncrement = Adduct.calculatePPMIncrement(monoisotopicMassPeak, monoisotopicMassAnnotation);
                    System.out.println("PPM Increment: " + ppmIncrement);
                    if (ppmIncrement <= tolerance) {
                        return entry.getKey();
                    }
                }
            }
        } else if (ionizationMode == IoniationMode.NEGATIVE) {
            for (Peak peak : groupedSignals) {
                double monoisotopicMassPeak = Adduct.getMonoisotopicMassFromMZ(peak.getMz(), "");
                for (Map.Entry<String, Double> entry : AdductList.MAPMZNEGATIVEADDUCTS.entrySet()) {
                    double monoisotopicMassAdduct = Adduct.getMonoisotopicMassFromMZ(mz, entry.getKey());
                    int ppmIncrement = Adduct.calculatePPMIncrement(monoisotopicMassPeak, monoisotopicMassAdduct);
                    if (ppmIncrement <= tolerance) {
                        return entry.getKey();
                    }
                }
            }
        }
        return "Unkown";
    }



    public static void main(String[] args) {

        // Given two peaks with ~21.98 Da difference (e.g., [M+H]+ and [M+Na]+)
        Peak mH = new Peak(700.500, 100000.0); // [M+H]+
        Peak mNa = new Peak(722.482, 80000.0);  // [M+Na]+
        Peak doublyCharged = new Peak(350.754, 85000.0);   // [M+2H]2+
        Lipid lipid = new Lipid(1, "PC 34:1", "C42H82NO8P", "PC", 34, 1);

        double annotationMZ = 700.49999d;
        double annotationMZNa = 722.47999d;
        double annotationMZdc = 350.75444d;
        double annotationIntensity = 80000.0;
        double annotationRT = 6.5d;
        Annotation annotation = new Annotation(lipid, annotationMZ, annotationIntensity, annotationRT, IoniationMode.POSITIVE, Set.of(mH, mNa, doublyCharged));
        Annotation annotation2 = new Annotation(lipid, annotationMZNa, annotationIntensity, annotationRT, IoniationMode.POSITIVE, Set.of(mH, mNa , doublyCharged));
        Annotation annotation3 = new Annotation(lipid, annotationMZdc, annotationIntensity, annotationRT, IoniationMode.POSITIVE, Set.of(mH, mNa, doublyCharged));



        // Then we should call the algorithmic/knowledge system rules fired to detect the adduct and Set it!
        //
        String adduct = annotation.detectAdductFromSignals(annotation.getIonizationMode());
        System.out.println("Adduct detected: " + adduct);
        annotation.setAdduct(adduct);
        System.out.println("Annotation: " + annotation);
        String adduct2 = annotation2.detectAdductFromSignals(annotation2.getIonizationMode());
        System.out.println("Adduct detected: " + adduct2);
        annotation2.setAdduct(adduct2);
        System.out.println("Annotation: " + annotation2);
        String adduct3 = annotation3.detectAdductFromSignals(annotation3.getIonizationMode());
        System.out.println("Adduct detected: " + adduct3);
        annotation3.setAdduct(adduct3);
        System.out.println("Annotation: " + annotation3);


    }





}






