import java.io.*;
import java.lang.invoke.SwitchPoint;

public class Main {

    private static String[] Motifs = new String[10];
    private static int Score = 0;
    private static int Iteration = 0;
    private static String deletedMotif = "";
    private static int deletedMotifIndex;
    private static String ConsensusString = "burda bi sorun var";

    public static void main(String[] args) {
//        CreateFile();
        String[] DNA = ReadFile();
        for (int m = 0; m < 10; m++) {
            for (int k = 9; k < 12; k++) {

                RandomizedMotifSearch(DNA, k);
                if (k == 9)
                    System.out.println("-----------\n|  k = " + k + "  |\n-----------");
                else
                    System.out.println("------------\n|  k = " + k + "  |\n------------");
                System.out.println("\nRANDOMIZED MOTIF SEARCH\n-----------------------\n\nMotifs :");
                for (int i = 0; i < 10; i++) {
                    System.out.println(Motifs[i]);
                }
                System.out.println("\nScore : " + Score + "\nConsensus String : " + ConsensusString);
                GibbsSampler(DNA, k);
                System.out.println("\nGIBBS SAMPLER\n-------------\n\nMotifs :");
                for (int i = 0; i < 10; i++) {
                    System.out.println(Motifs[i]);
                }
                System.out.println("\nScore : " + Score + "\nConsensus String : " + ConsensusString);
                System.out.println("-----------------------------------------------------------------------");
            }
        }
    }

    /*    READ WRITE    */
    public static void CreateFile() {

        String str = "";
        for (int i = 0; i < 10; i++) {
            str += createRandomString();
            str += "\n";
        }
        String tenMer = createTenMer();
        String[] allTenMers = new String[10];
        System.out.println("Normal 10-mer : \n" + tenMer);
        System.out.println("All 10-mers : ");
        for (int i = 0; i < 10; i++) {
            String mutationStr = tenMer;
            int[] indexs = MutationIndexs();
            for (int j = 0; j < indexs.length; j++) {
                mutationStr = Mutation(mutationStr, indexs[j]);
            }
            allTenMers[i] = mutationStr;
        }
        for (int i = 0; i < 10; i++) {
            int index = (int) (Math.random() * 490);
            System.out.print(allTenMers[i] + " ");
            System.out.println(index);
            index = i * 500 + index;
            str = str.substring(0, index) + allTenMers[i] + str.substring(index + 10);
        }

        File file = new File("C:/Users/pelsi/Desktop/input.txt");
        if (!file.exists())
            try {
                file.createNewFile();
            } catch (Exception e) {
                e.printStackTrace();
            }

        try {
            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write(str);
            bw.close();
            System.out.println("SUCCESS..");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static String[] ReadFile() {
        String[] motifs = new String[10];
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(
                    "C:/Users/pelsi/Desktop/input.txt"));
            String line = reader.readLine();
            int i = 0;
            while (line != null) {
                motifs[i] = line;
                i++;
                // read next line
                line = reader.readLine();
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return motifs;
    }

    /*    ALGORITHMS   */

    public static void RandomizedMotifSearch(String[] DNA, int k) {
        String[] motifs = selectRandomMotifs(DNA, k);
        String[] bestMotifs = motifs.clone();
        int iteration = 0, counter = 0;
        while (true) {
            iteration++;
            double[][] profile = calculateProfile(motifs, k);
            motifs = selectNewMotifs(DNA, profile, k);
            int motifScore = calculateScore(motifs, k);
            int bestMotifScore = calculateScore(bestMotifs, k);
            if (motifScore <= bestMotifScore) {
                bestMotifs = motifs.clone();
                if (motifScore == bestMotifScore)
                    counter++;
                if (counter > 50) {
                    Motifs = bestMotifs;
                    Score = bestMotifScore;
                    Iteration = iteration;
                    ConsensusString = getConsensusString(Motifs, profile, k);
                    break;
                }
            } else {
                Motifs = bestMotifs;
                Score = bestMotifScore;
                Iteration = iteration;
                ConsensusString = getConsensusString(Motifs, profile, k);
                break;
            }
        }

    }

    public static void GibbsSampler(String[] DNA, int k) {
        String[] motifs = selectRandomMotifs(DNA, k);
        int iteration = 0;
        int counter = 0;
        while (true) {
            iteration++;
            int scoreMotifs = calculateScore(motifs, k);
            String[] newMotifsDeletedOne = deleteRandomMotif(motifs);
            double[][] profile = calculateProfileGibbs(newMotifsDeletedOne, k);
            String newMotif = selectMotifFromDeletedLine(DNA, k, profile);
            String[] newMotifs = putNewMotifToMotifs(newMotifsDeletedOne, newMotif);
            int scoreNewMotifs = calculateScore(newMotifs, k);
            if (scoreMotifs >= scoreNewMotifs) {
                motifs = newMotifs;
                if (scoreMotifs == scoreNewMotifs)
                    counter++;
                if (counter > 50) {
                    Iteration = iteration;
                    Score = scoreNewMotifs;
                    Motifs = newMotifs;
                    ConsensusString = getConsensusString(Motifs, profile, k);
                    break;
                }
            } else {
                Iteration = iteration;
                Score = scoreNewMotifs;
                Motifs = newMotifs;
                ConsensusString = getConsensusString(Motifs, profile, k);
                break;
            }

        }

    }

    /*    OTHER METHODS   */

    public static String[] putNewMotifToMotifs(String[] motifs, String newMotif) {
        String[] newMotifs = new String[10];
        for (int i = 0; i < 10; i++) {
            if (i < deletedMotifIndex) {
                newMotifs[i] = motifs[i];
            } else if (deletedMotifIndex < i) {
                newMotifs[i] = motifs[i - 1];
            } else {
                newMotifs[i] = newMotif;
            }
        }
        return newMotifs;
    }

    public static String selectMotifFromDeletedLine(String[] DNA, int k, double[][] profile) {
        String deletedDNA = DNA[deletedMotifIndex];
        double maxProb = 0;
        int index = 0;
        for (int i = 0; i < 500 - k; i++) {
            String motif = deletedDNA.substring(i, i + k);
            double prob = calculateProfileScore(motif, profile);
            if (prob > maxProb) {
                maxProb = prob;
                index = i;
            }
        }
        String newMotif = deletedDNA.substring(index, index + k);
        return newMotif;
    }

    public static String[] deleteRandomMotif(String[] motifs) {
        String[] newMotifs = new String[9];
        int a = (int) (Math.random() * 10);
        deletedMotif = motifs[a];
        deletedMotifIndex = a;
        for (int i = 0, k = 0; i < motifs.length; i++) {
            if (i == a) {
                continue;
            }
            newMotifs[k++] = motifs[i];
        }
        return newMotifs;
    }

    public static int calculateScore(String[] motifs, int k) {
        int score = 0;
        for (int i = 0; i < k; i++) {
            int a = 0, c = 0, g = 0, t = 0;
            for (int j = 0; j < 10; j++) {
                if (motifs[j].charAt(i) == 'a')
                    a++;
                else if (motifs[j].charAt(i) == 'c')
                    c++;
                else if (motifs[j].charAt(i) == 'g')
                    g++;
                else if (motifs[j].charAt(i) == 't')
                    t++;
            }
            int max = Math.max(Math.max(a, c), Math.max(g, t));
            score += 10 - max;
        }
        return score;
    }

    public static String[] selectNewMotifs(String[] DNA, double[][] profile, int k) {
        String[] newMotifs = new String[10];
        for (int i = 0; i < 10; i++) {
            double max = 0;
            String motif = "";
            for (int j = 0; j < 500 - k; j++) {
                double calculatedProfile = calculateProfileScore(DNA[i].substring(j, j + k), profile);
                if (calculatedProfile > max) {
                    max = calculatedProfile;
                    motif = DNA[i].substring(j, j + k);
                }
            }
            newMotifs[i] = motif;
        }
        return newMotifs;
    }

    public static double calculateProfileScore(String motif, double[][] profile) {
        double profileScore = 1;
        for (int i = 0; i < motif.length(); i++) {
            if (motif.charAt(i) == 'a')
                profileScore *= profile[0][i];
            else if (motif.charAt(i) == 'c')
                profileScore *= profile[1][i];
            else if (motif.charAt(i) == 'g')
                profileScore *= profile[2][i];
            else if (motif.charAt(i) == 't')
                profileScore *= profile[3][i];
            if (profileScore == 0)
                break;
        }
        return profileScore;
    }

    public static double[][] calculateProfile(String[] motifs, int k) {
        double[][] profile = new double[4][k];
        for (int i = 0; i < k; i++) {
            int a = 0, c = 0, g = 0, t = 0;
            for (String motif : motifs) {
                if (motif.charAt(i) == 'a')
                    a++;
                else if (motif.charAt(i) == 'c')
                    c++;
                else if (motif.charAt(i) == 'g')
                    g++;
                else if (motif.charAt(i) == 't')
                    t++;
            }
            profile[0][i] = (double) a / 10;
            profile[1][i] = (double) c / 10;
            profile[2][i] = (double) g / 10;
            profile[3][i] = (double) t / 10;
        }
        return profile;
    }

    public static double[][] calculateProfileGibbs(String[] motifs, int k) {
        double[][] profile = new double[4][k];
        for (int i = 0; i < k; i++) {
            int a = 1, c = 1, g = 1, t = 1;
            for (String motif : motifs) {
                if (motif.charAt(i) == 'a')
                    a++;
                else if (motif.charAt(i) == 'c')
                    c++;
                else if (motif.charAt(i) == 'g')
                    g++;
                else if (motif.charAt(i) == 't')
                    t++;
            }
            double divide = motifs.length + 4;
            profile[0][i] = (double) a / divide;
            profile[1][i] = (double) c / divide;
            profile[2][i] = (double) g / divide;
            profile[3][i] = (double) t / divide;
        }
        return profile;
    }

    public static String[] selectRandomMotifs(String[] DNA, int k) {
        String[] motifs = new String[10];
        for (int i = 0; i < 10; i++) {
            int random = (int) (Math.random() * (500 - k));
            motifs[i] = DNA[i].substring(random, random + k);
        }
        return motifs;
    }

    public static String generateRandomNucleotide() {
        String nucleotid = "";
        int a = (int) (Math.random() * 4);
        switch (a) {
            case 0:
                nucleotid = "a";
                break;
            case 1:
                nucleotid = "g";
                break;
            case 2:
                nucleotid = "t";
                break;
            case 3:
                nucleotid = "c";
                break;
            default:
        }
        return nucleotid;
    }

    public static String createTenMer() {
        String str = "";
        for (int i = 0; i < 10; i++) {
            str = str + generateRandomNucleotide();
        }
        return str;
    }

    public static String createRandomString() {
        String str = "";
        for (int i = 0; i < 500; i++) {
            str = str + generateRandomNucleotide();
        }
        return str;

    }

    public static String Mutation(String str, int index) {
        while (true) {
            String nucleotid = generateRandomNucleotide();
            if (nucleotid.charAt(0) != str.charAt(index)) {
                str = str.substring(0, index) + nucleotid + str.substring(index + 1);
                break;
            }
        }
        return str;
    }

    public static int[] MutationIndexs() {
        int[] indexs = new int[4];
        int a = (int) (Math.random() * 10);
        indexs[0] = a;
        for (int i = 1; i < 4; i++) {
            a = (int) (Math.random() * 10);
            if (indexs[i - 1] != a) {
                indexs[i] = a;
            } else
                i--;
        }
        return indexs;
    }

    public static String getConsensusString(String[] motifs, double[][] profile, int k) {
        String consensusString = "";
        for (int i = 0; i < k; i++) {
            int a = 0, c = 0, g = 0, t = 0;
            for (int j = 0; j < 10; j++) {
                if (motifs[j].charAt(i) == 'a')
                    a++;
                else if (motifs[j].charAt(i) == 'c')
                    c++;
                else if (motifs[j].charAt(i) == 'g')
                    g++;
                else if (motifs[j].charAt(i) == 't')
                    t++;
            }
            int max = Math.max(Math.max(a, c), Math.max(g, t));
            if (a == max)
                consensusString += "a";
            else if (c == max)
                consensusString += "c";
            else if (g == max)
                consensusString += "g";
            else if (t == max)
                consensusString += "t";
            else
                consensusString += "?";
        }
        return consensusString;
    }
}
