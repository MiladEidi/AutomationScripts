/**
 * @author Milad
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class VCFsubsetter {

    public static void main(String[] args) {
        String inputFilePath = "G:\\Eidi\\ExAC.r1.sites.vep.vcf"; // Path to your input VCF file
        int variantsPerSubset = 700000; // Number of variants per subset VCF file

        try {
            subsetVCF(inputFilePath, variantsPerSubset);
            System.out.println("VCF file subsetted successfully.");
        } catch (IOException e) {
            System.err.println("Error processing VCF file: " + e.getMessage());
        }
    }

    public static void subsetVCF(String inputFilePath, int variantsPerSubset) throws IOException {
        List<String> headerLines = new ArrayList<>();

        // Read the header lines
        try (BufferedReader reader = new BufferedReader(new FileReader(inputFilePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) { // Header line
                    headerLines.add(line);
                } else {
                    break; // Stop reading after reaching the end of header
                }
            }
        }

        // Create subset VCF files
        try (BufferedReader reader = new BufferedReader(new FileReader(inputFilePath))) {
            String outputFilePath = null;
            BufferedWriter writer = null;
            int variantsWritten = 0;
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    // Skip header lines as they are already saved
                    continue;
                }
                if (variantsWritten % variantsPerSubset == 0) {
                    // Close previous writer and open new file for next subset
                    if (writer != null) {
                        writer.close();
                    }
                    outputFilePath = "G:\\Eidi\\subset_" + ((variantsWritten / variantsPerSubset) + 1) + ".vcf";
                    writer = new BufferedWriter(new FileWriter(outputFilePath));
                    // Write header lines
                    for (String headerLine : headerLines) {
                        writer.write(headerLine);
                        writer.newLine();
                    }
                }
                // Write variant line
                writer.write(line);
                writer.newLine();
                variantsWritten++;
            }
            // Close the last writer
            if (writer != null) {
                writer.close();
            }
        }
    }
}

