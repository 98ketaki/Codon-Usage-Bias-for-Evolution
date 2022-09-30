package main

import (
	"bufio"
	"fmt"
	"image/color"
	"log"
	"os"
	"strings"

	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg/draw"
)

type xy struct {
	x []float64
	y []float64
}

//OpenInputFile opens the input file
//It takes the fileName as a string and returns the lines in the file as a list of strings
func OpenInputFile(fileName string) []string {
	seq, err := os.Open(fileName)

	if err != nil {
		fmt.Println("Error: Cannot open the input file")
		os.Exit(3)
	}

	// Create the variable to hold the lines
	var lines []string = make([]string, 0)

	scanner := bufio.NewScanner(seq)

	for scanner.Scan() {
		// Append it to the lines slice
		lines = append(lines, scanner.Text())
	}

	if scanner.Err() != nil {
		fmt.Println("Sorry: Error during the file reading")
		os.Exit(3)
	}

	// close the file and return the lines
	seq.Close()
	return lines
}

//GetMapAccessionID gives the Accession ID of the sequences
//It takes the data from file as a list of strings and returns the accession IDs as a list of strings
func GetMapAccessionID(fileData []string) []string {
	var accessID []string = make([]string, 0)

	for i := 0; i < len(fileData); i++ {
		startTag := strings.Contains(fileData[i], ">")

		if startTag == true {
			var data []string = strings.Split(fileData[i], " ")
			accessionID := data[0]
			accessionID = accessionID[1:len(accessionID)]
			accessID = append(accessID, accessionID)
		}
	}

	return accessID
}

//GetNucleotideSequence gives the sequences from the input file
//It takes the data from file as a list of strings and returns the individual sequences as a list of strings
func GetNucleotideSequence(fileData []string) []string {

	var nucleotideSequence []string = make([]string, 0)
	var singleNucleotideSequence string

	for i := 0; i < len(fileData); i++ {
		StartTag := strings.Contains(fileData[i], ">")

		if StartTag == false {

			if fileData[i] == "" {
				nucleotideSequence = append(nucleotideSequence, singleNucleotideSequence)
				singleNucleotideSequence = ""
			} else {
				singleNucleotideSequence = singleNucleotideSequence + fileData[i]
			}
		}
	}

	nucleotideSequence = append(nucleotideSequence, singleNucleotideSequence)
	singleNucleotideSequence = ""

	return nucleotideSequence
}

//Mapper maps the accession IDs to their nucleotide sequences
//It takes a list of Accession ID and list of sequences and maps each ID to its sequence
func Mapper(listOfAccessionID, nucleotideSequence []string) map[string]string {

	mappedIdToNucleotideSeq := make(map[string]string)

	accessionIDLen := len(listOfAccessionID)
	nucleotideSequenceLen := len(nucleotideSequence)

	if accessionIDLen != nucleotideSequenceLen {
		fmt.Println("Error: the number of accession IDs and sequences are not matching")
		os.Exit(3)
	} else {
		for i := 0; i < nucleotideSequenceLen; i++ {

			mappedIdToNucleotideSeq[listOfAccessionID[i]] = nucleotideSequence[i]
		}
	}

	return mappedIdToNucleotideSeq
}

//CreateFile creates a file
//It creates the output file
func CreateFile(fileName string) {
	// Open the file to write the frequency map
	outFile, err := os.Create(fileName)

	if err != nil {
		fmt.Println("Sorry: Could not create the file")
	}

	outFile.Close()
}

//WriteDataToFile writes results in the output file
//It writes the RSCU, GC, Neutrality and Parity Data in the output file
func WriteDataToFile(fileName, accessionID string, freqMap map[string]float64, xAxisNeutralityplot, yAxisNeutralityplot, xAxisParityPlot, yAxisParityPlot, totalNuclearFrequencyfloat float64) {
	outFile, err := os.OpenFile(fileName, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		fmt.Println("Error: Cannot open the input files")
		os.Exit(3)
	}

	fmt.Fprintf(outFile, "RSCU data for accession ID : %v \n", accessionID)

	for sequence := range freqMap {
		fmt.Fprintf(outFile, "%v : %f \n", sequence, freqMap[sequence])
	}
	fmt.Fprintf(outFile, "X Axis for neutralityplot %f,Y Axis for neutralityplot %f  \n", xAxisNeutralityplot, yAxisNeutralityplot)
	fmt.Fprintf(outFile, "X Axis for parityplot %f,Y Axis for parityplot %f  \n", xAxisParityPlot, yAxisParityPlot)
	fmt.Fprintf(outFile, "GC content is  %f \n", totalNuclearFrequencyfloat)
	fmt.Fprintf(outFile, "\n \n")
}

//WriteDataPoints writes all the X and Y co-ordinates for neutrality and Parity analyses to their respective text files
// It takes the x and y co-ordinate values for each data points and the name of output file as a string
func WriteDataPoints(x, y float64, fileName string) {
	outFile, err := os.OpenFile(fileName, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)

	if err != nil {
		fmt.Println("Error: Cannot open the output file")
		os.Exit(3)
	}

	fmt.Fprintf(outFile, "%f,%f \n", x, y)

}

//ReadData takes the outputfile name and returns the x and y co-ordinates
//It opens the neutrality and parity output files and returns the  x and y co-ordinates for plotting
func ReadData(path string) (plotter.XYs, float64, float64, error) {
	f, err := os.Open(path)

	if err != nil {
		return nil, 0, 0, err
	}
	defer f.Close()

	var listofx []float64
	var listofy []float64

	data := xy{

		x: listofx,

		y: listofy,
	}

	var xys plotter.XYs
	s := bufio.NewScanner(f)

	for s.Scan() {
		var x, y float64
		_, err := fmt.Sscanf(s.Text(), "%f,%f", &x, &y)
		if err != nil {
			log.Printf("Discarding bad data point %q: %v", s.Text(), err)
			continue
		}
		xys = append(xys, struct{ X, Y float64 }{x, y})
		data.y = append(data.y, y)
		data.x = append(data.x, x)
	}

	if err := s.Err(); err != nil {
		return nil, 0, 0, fmt.Errorf("Error: Could not scan the x and y values: %v", err)
	}

	//Regression Line for the Neutrality Plot
	b, a := stat.LinearRegression(data.x, data.y, nil, false)

	return xys, b, a, nil
}

//plotData takes the output filename and creates the plot with that name
//It takes name of the output file and the co-ordinates
//It uses the gonum plot package to plot the scatter plot
//The regression line is plotted
//It also determines the mutational pressure from the data and prints it
func plotData(path string, xys plotter.XYs, b, a float64, flag int) error {
	f, err := os.Create(path)

	if err != nil {
		return fmt.Errorf("Cannot create the output file for the plot%s: %v", path, err)
	}

	p := plot.New()
	if err != nil {
		return fmt.Errorf("Cannot create plot: %v", err)
	}

	// create scatter with all data points
	s, err := plotter.NewScatter(xys)

	if err != nil {
		return fmt.Errorf("Cannot create scatter plot: %v", err)
	}
	s.GlyphStyle.Shape = draw.CrossGlyph{}
	s.Color = color.RGBA{R: 255, A: 255}
	p.Add(s)

	//Neutrality Plot
	if flag == 1 {
		fmt.Println("The effect of mutational pressure is", a)
		// y = ax + b ; have a,b
		// create a linear regression line using a and b values obtained
		l, err := plotter.NewLine(plotter.XYs{
			{0.3, 0.3*a + b}, {0.9, 0.9*a + b},
		})

		p.Title.Text = "Neutrality Plot"
		p.X.Label.Text = "GC3"
		p.Y.Label.Text = "GC12"

		if err != nil {
			return fmt.Errorf("Cannot create the regression line: %v", err)
		}
		p.Add(l)
	}

	//Parity Plot
	if flag == 0 {
		p.Title.Text = "Parity Plot"
		p.X.Label.Text = "G3/GC3"
		p.Y.Label.Text = "A3/AT3"
	}

	wt, err := p.WriterTo(256, 256, "png")
	if err != nil {
		return fmt.Errorf("Cannot create writer: %v", err)
	}

	_, err = wt.WriteTo(f)
	if err != nil {
		return fmt.Errorf("Cannot write to %s: %v", path, err)
	}

	if err := f.Close(); err != nil {
		return fmt.Errorf("Cannot close the plot %s: %v", path, err)
	}

	return nil
}

//CodonFinder returns the frequency of each codon in the sequences
// It takes the sequence a input and returns the frequency of each codon
func CodonFinder(seq string) map[string]int {
	length := len(seq)
	var a []string

	for i := 0; i < length-2; i += 3 {
		codon := (seq[i : i+3])
		a = append(a, codon)
	}

	Frequency := CodonFrequency(a)
	return Frequency
}

//CodonFrequency calculates frequency of each codon
//It takes a list of codons from the sequence to return the frequency
func CodonFrequency(codons []string) map[string]int {
	freq := make(map[string]int)

	for _, i := range codons {
		if freq[i] == 0 {
			freq[i] = 1
		} else {
			freq[i]++
		}
	}

	return freq
}

//RscuCalc takes the codon frequency and list of synonymous codons that code each amino acid to
//return the relative synonymous codon usage (RSCU) values

//Relative Synonymous Codon Usage value for a Codon can be calculated using
//the observed number of occurrences  divided by the number expected if all synonymous codons
//for an amino acid were used equally frequently
func RscuCalc(codonFreq map[string]int, synonymousCodons map[string][]string) map[string]float64 {
	rscu := make(map[string]float64)

	for aminoAcid := range synonymousCodons {
		total := 0
		codons := synonymousCodons[aminoAcid]

		for _, codon := range codons {
			total += codonFreq[codon]
		}

		if total != 0 {
			for _, codon := range codons {
				d := float64(total) / float64(len(codons))
				rscuValue := float64(codonFreq[codon]) / d
				rscu[codon] = rscuValue
			}
		}

	}

	return rscu
}

//TotalNucFreq takes the input sequence to return the GC content
//GC content gives the proportion of the nucleotides G and C in the given sequence
func TotalNucFreq(seq string) float64 {
	a := len(seq)
	var seq1 []string

	for x := 0; x < a; x = x + 1 {
		z := string(seq[x])
		seq1 = append(seq1, z)
	}
	freq := make(map[string]int)

	for _, i := range seq1 {
		if freq[i] == 0 {
			freq[i] = 1
		} else {
			freq[i]++
		}
	}

	return GcContent(freq)
}

//PartialNucContent takes the input sequence and returns the nucleotide content at each
//position of the codon
func PartialNucContent(seq string) (map[string]int, map[string]int, map[string]int) {
	a := len(seq)
	var seq1 []string
	var seq2 []string
	var seq3 []string

	for x := 0; x < a; x = x + 3 {
		z := string(seq[x])
		seq1 = append(seq1, z)
	}

	for x := 1; x < a; x = x + 3 {
		z := string(seq[x])
		seq2 = append(seq2, z)
	}

	for x := 2; x < a; x = x + 3 {
		z := string(seq[x])
		seq3 = append(seq3, z)
	}

	return PartialNucFreq(seq1, seq2, seq3)
}

//PartialNucFreq takes the nucleotide at each codon position as input and
//return the nucleotide frequency at each codon position
func PartialNucFreq(seq1, seq2, seq3 []string) (map[string]int, map[string]int, map[string]int) {
	freq1 := make(map[string]int)
	freq2 := make(map[string]int)
	freq3 := make(map[string]int)

	for _, i := range seq1 {
		if freq1[i] == 0 {
			freq1[i] = 1
		} else {
			freq1[i]++
		}
	}

	for _, j := range seq2 {
		if freq2[j] == 0 {
			freq2[j] = 1
		} else {
			freq2[j]++
		}
	}

	for _, k := range seq3 {
		if freq3[k] == 0 {
			freq3[k] = 1
		} else {
			freq3[k]++
		}
	}
	//n := NeutralityPlot(freq1, freq2, freq3)
	return freq1, freq2, freq3
	//return ParityPlot(freq3)
}

//GcContent takes the frequency of nucleotide and returns the GC content of the sequence
func GcContent(gcFreq map[string]int) float64 {
	var a = float64(gcFreq["A"])
	var c = float64(gcFreq["C"])
	var g = float64(gcFreq["G"])
	var t = float64(gcFreq["T"])
	gcContent := ((g + c) / (a + c + g + t))
	return gcContent
}

//PartialGcContent takes the frequency of nucleotide at each position and returns the GC content of the sequence
func PartialGcContent(gcFreq map[string]int) float64 {
	var a = float64(gcFreq["A"])
	var c = float64(gcFreq["C"])
	var g = float64(gcFreq["G"])
	var t = float64(gcFreq["T"])
	gcContent := ((g + c) / (a + c + g + t))
	return gcContent
}

//NeutralityPlot takes the frequency of nucleotides at each codon position and returns the values of gc12 and gc3
//gc12 is the avaerage gc content at first and second codon position
//gc3 is the gc content at third codon position
//Neutrality plot is a scatter plot between gc12 on Y axis and gc3 on X axis
func NeutralityPlot(freq1, freq2, freq3 map[string]int) (float64, float64) {
	gc1 := (float64(freq1["G"]) + float64(freq1["C"])) / (float64(freq3["G"]) + float64(freq3["C"]) + float64(freq3["A"]) + float64(freq3["T"]))
	gc2 := (float64(freq2["G"]) + float64(freq2["C"])) / (float64(freq3["G"]) + float64(freq3["C"]) + float64(freq3["A"]) + float64(freq3["T"]))
	gc12 := (gc1 + gc2) / float64(2)
	gc3 := (float64(freq3["G"]) + float64(freq3["C"])) / (float64(freq3["G"]) + float64(freq3["C"]) + float64(freq3["A"]) + float64(freq3["T"]))
	return gc12, gc3
}

//ParityPlot takes the frequency of nucleotides at third codon position and returns the values of a3, g3, at3 and gc3
//a3 and g3 is the proportion of adenine and guanine at third codon position
//at3 and gc3 is the proportion of adenine and thymine,  guanine and cytosine respectively at third codon position
//Parity plot is a scatter plot between a3 / at3 on Y axis and g3 / gc3 on X axis
func ParityPlot(freq3 map[string]int) (float64, float64) {
	a3 := float64(freq3["A"])
	g3 := float64(freq3["G"])
	gc3 := float64(freq3["G"]) + float64(freq3["C"])
	at3 := float64(freq3["A"]) + float64(freq3["T"])
	xAxis := a3 / at3
	yAxis := g3 / gc3
	return xAxis, yAxis
}

//SynonymousCodons has a map of each amino acid and its corresponding codons
func SynonymousCodons() map[string][]string {
	synonymousCodons := make(map[string][]string)
	synonymousCodons["CYS"] = append(synonymousCodons["CYS"], "TGT", "TGC")
	synonymousCodons["ASP"] = append(synonymousCodons["ASP"], "GAT", "GAC")
	synonymousCodons["SER"] = append(synonymousCodons["SER"], "TCT", "TCG", "TCA", "TCC", "AGC", "AGT")
	synonymousCodons["GLN"] = append(synonymousCodons["GLN"], "CAA", "CAG")
	synonymousCodons["MET"] = append(synonymousCodons["MET"], "ATG")
	synonymousCodons["ASN"] = append(synonymousCodons["ASN"], "AAC", "AAT")
	synonymousCodons["PRO"] = append(synonymousCodons["PRO"], "CCT", "CCG", "CCA", "CCC")
	synonymousCodons["LYS"] = append(synonymousCodons["LYS"], "AAG", "AAA")
	synonymousCodons["THR"] = append(synonymousCodons["THR"], "ACC", "ACA", "ACG", "ACT")
	synonymousCodons["PHE"] = append(synonymousCodons["PHE"], "TTT", "TTC")
	synonymousCodons["ALA"] = append(synonymousCodons["ALA"], "GCA", "GCC", "GCG", "GCT")
	synonymousCodons["GLY"] = append(synonymousCodons["GLY"], "GGT", "GGG", "GGA", "GGC")
	synonymousCodons["ILE"] = append(synonymousCodons["ILE"], "ATC", "ATA", "ATT")
	synonymousCodons["LEU"] = append(synonymousCodons["LEU"], "TTA", "TTG", "CTC", "CTT", "CTG", "CTA")
	synonymousCodons["HIS"] = append(synonymousCodons["HIS"], "CAT", "CAC")
	synonymousCodons["ARG"] = append(synonymousCodons["ARG"], "CGA", "CGC", "CGG", "CGT", "AGG", "AGA")
	synonymousCodons["TRP"] = append(synonymousCodons["TRP"], "TGG")
	synonymousCodons["VAL"] = append(synonymousCodons["VAL"], "GTA", "GTC", "GTG", "GTT")
	synonymousCodons["GLU"] = append(synonymousCodons["GLU"], "GAG", "GAA")
	synonymousCodons["TYR"] = append(synonymousCodons["TYR"], "TAT", "TAC")
	return synonymousCodons
}

func main() {

	//Access the Codon Data and data handling functions
	filename := os.Args[1]                                          //Read the user given filename
	fileData := OpenInputFile(filename)                             //Open the file with given filename and return the file Data
	listOfAccessionID := GetMapAccessionID(fileData)                //Get accession ID from data file
	nucleotideSequence := GetNucleotideSequence(fileData)           //Get the nucleotide sequence
	mappedSequence := Mapper(listOfAccessionID, nucleotideSequence) //Map accession ID to the nucleotide sequence

	//Create synonymous codons list
	synonymousCodons := SynonymousCodons()

	//Create the required output file and create the input data file for the Parity graph and neutrality graph
	fileName := "output.txt"
	parityFile := "Parity.txt"
	neutralityFile := "Neutrality.txt"

	CreateFile(fileName)
	CreateFile(parityFile)
	CreateFile(neutralityFile)

	//Calculate various values based on the input NucleotideSequence against the synonymous codon data
	for i := 0; i < len(listOfAccessionID); i++ {

		nucleotideSequence := mappedSequence[listOfAccessionID[i]]

		codonFreq := CodonFinder(nucleotideSequence)
		rscuMappedData := RscuCalc(codonFreq, synonymousCodons)

		freq1, freq2, freq3 := PartialNucContent(nucleotideSequence)

		xAxisNeutralityplot, yAxisNeutralityplot := NeutralityPlot(freq1, freq2, freq3)

		xAxisParityPlot, yAxisParityPlot := ParityPlot(freq3)

		TotalNuclearFrequency := TotalNucFreq(nucleotideSequence)

		//Write all the resulting Value to their respective text files
		WriteDataPoints(xAxisNeutralityplot, yAxisNeutralityplot, neutralityFile)
		WriteDataPoints(xAxisParityPlot, yAxisParityPlot, parityFile)
		WriteDataToFile(fileName, listOfAccessionID[i], rscuMappedData, xAxisNeutralityplot, yAxisNeutralityplot, xAxisParityPlot, yAxisParityPlot, TotalNuclearFrequency)
	}

	//Create the datapoints structure  and create the graph
	//Flag value depicts if regression line is present on final graph or not
	//b,a is data returned from the function which depicts regression values in the equation y = ax + b
	//Creating the Neutrality Plot
	xys, b, a, err := ReadData("Neutrality.txt")

	if err != nil {
		log.Fatalf("could not read data.txt: %v", err)
	}
	_ = xys
	flag := 1
	err = plotData("NeutralityGraph.png", xys, b, a, flag) // Regression Line is present

	if err != nil {
		log.Fatalf("could not plot data: %v", err)
	}

	//Creating the Parity Plot
	xys, b, a, err = ReadData("Parity.txt")

	if err != nil {
		log.Fatalf("could not read data.txt: %v", err)
	}
	_ = xys

	flag = 0
	err = plotData("ParityGraph.png", xys, b, a, flag) // Regression Line is not present

	if err != nil {
		log.Fatalf("could not plot data: %v", err)
	}

}
