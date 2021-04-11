import java.io.*;

public class Main {

    private static int[] indexOfKMers = new int[10];

    public static void main(String[] args) {
        //WriteFile();
        String[] DNA = ReadFile();
        for (int i = 0; i < DNA.length; i++) {
            System.out.println(DNA[i]);
        }
        RandomizedMotifSearch(DNA, 10);
    }

    public static void WriteFile() {

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

        //String oldStr = str;
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


    public static void RandomizedMotifSearch(String[] DNA, int k) {
        String[] motifStrings = new String[10];
        int[] indexs = new int[10];
        int[] bestIndexs;
        for (int i = 0; i < 10; i++) {
            int a = (int) (Math.random() * (500 - k));
            motifStrings[i] = DNA[i].substring(a, a + k);
            indexs[i] = a;
        }
        String[] bestMotifStrings = motifStrings.clone();

        int[][] motifs;
        double[][] profile;
        int[][] bestMotifs;
        double[][] bestProfile;
        int counter = 0;
        while (true) {
            bestMotifs = Motifs(bestMotifStrings);
            bestProfile = Profile(bestMotifs);
            bestIndexs = UseProfileInDNA(bestProfile, DNA, k);

            for (int i = 0; i < 10; i++)
                motifStrings[i] = DNA[i].substring(indexs[i], indexs[i] + k);
            motifs = Motifs(motifStrings);
            profile = Profile(motifs);
            indexs = UseProfileInDNA(profile, DNA, k);

            if (ScoreAllMotifs(motifStrings, profile) < ScoreAllMotifs(bestMotifStrings, bestProfile)) {
                bestMotifStrings = motifStrings.clone();
                bestIndexs = indexs;
            } else
//                counter++;
//            if (counter > 50)
                break;
        }
        for (int i = 0; i < 10; i++) {
            System.out.println(bestMotifStrings[i] + " " + bestIndexs[i]);
        }
    }

    public static double Score(String motif, double[][] profile) {
        double score = 1;
        for (int i = 0; i < motif.length(); i++) {
            if (motif.charAt(i) == 'a') {
                score *= profile[0][i];
            } else if (motif.charAt(i) == 'c') {
                score *= profile[1][i];
            } else if (motif.charAt(i) == 'g') {
                score *= profile[2][i];
            } else if (motif.charAt(i) == 't') {
                score *= profile[3][i];
            }
            if (score == 0)
                return 0;
        }
        return score;
    }

    public static double ScoreAllMotifs(String[] motifs, double[][] profile) {
        double sumScore = 0;
        double score = 1;
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < motifs[j].length(); i++) {
                if (motifs[j].charAt(i) == 'a') {
                    score *= profile[0][i];
                } else if (motifs[j].charAt(i) == 'c') {
                    score *= profile[1][i];
                } else if (motifs[j].charAt(i) == 'g') {
                    score *= profile[2][i];
                } else if (motifs[j].charAt(i) == 't') {
                    score *= profile[3][i];
                }
                if (score == 0)
                    break;
            }
            sumScore += score;
        }
        return sumScore;
    }

    public static int[][] Motifs(String[] motifStrings) {
        int[][] motifs = new int[4][10];
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                if (motifStrings[i].charAt(j) == 'a') {
                    motifs[0][j]++;
                } else if (motifStrings[i].charAt(j) == 'c') {
                    motifs[1][j]++;
                } else if (motifStrings[i].charAt(j) == 'g') {
                    motifs[2][j]++;
                } else if (motifStrings[i].charAt(j) == 't') {
                    motifs[3][j]++;
                }
            }
        }
        return motifs;
    }

    public static double[][] Profile(int[][] motifs) {
        double[][] profile = new double[4][10];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 10; j++) {
                profile[i][j] = (double) motifs[i][j] / (double) 10;
            }
        }
        return profile;
    }

    public static int[] UseProfileInDNA(double[][] profile, String[] DNA, int k) {
        int[] indexs = new int[10];
        double score = 0;
        for (int i = 0; i < 10; i++) {
            int counter = 0;
            while (counter < 500 - k) {
                String controlKMer = DNA[i].substring(counter, counter + k);
                double newScore = Score(controlKMer, profile);
                if (newScore > score) {
                    score = newScore;
                    indexs[i] = counter;
                }
                counter++;
            }
        }
        return indexs;
    }


    public static String createTenMer() {
        String str = "";
        for (int i = 0; i < 10; i++) {
            str = str + generateRandomNucleotid();
        }
        return str;
    }

    public static String createRandomString() {
        String str = "";
        for (int i = 0; i < 500; i++) {
            str = str + generateRandomNucleotid();
        }
        return str;

    }

    public static String Mutation(String str, int index) {
        while (true) {
            String nucleotid = generateRandomNucleotid();
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

    public static String generateRandomNucleotid() {
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

}
