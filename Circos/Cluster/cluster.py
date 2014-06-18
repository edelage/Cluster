"""
This program creates files needed by Circos in order to link Differentially Expressed (DE) genes
found by DESeq2 with clusters given by WGCNA analysis.
"""

import sys
import argparse

## Global variables
contenuCluster = "" ## Stores gene-module association
contenuCaryotype = "" ## Stores info about the modules
contenuGeneDiff = "" ## Stores info about the differentially expressed genes

def main():
    
    argParser = argparse.ArgumentParser()
    argParser.add_argument("geneModule",help="List the gene-module association given by WGCNA")
    argParser.add_argument("geneDiff",help="List of differentially expressed genes")
    argParser.add_argument("corrModule",help="Give the correlation between a module and a specific trait")
    argParser.add_argument("dirData",help="Directory where temp files are generated")


    # If no arguments were given, print the help
    if len(sys.argv) == 1:
        argParser.print_help()
        sys.exit(1)
        
    args = argParser.parse_args()
    
    ## Constant definition
    outFileCluster = args.dirData + "/cluster.txt"
    outFileCaryotype = args.dirData + "/caryotype.txt"
    outFileGeneDiff = args.dirData + "/gene.names.txt" 
    outFileHeatmap = args.dirData + "/heatmapage.txt" 
    outFileHistogram = args.dirData + "/histo"                   
   

    ## Program core
    buildClusterCaryotype(args.geneModule,outFileCluster,outFileCaryotype)
    buildGeneDiff(args.geneDiff,outFileGeneDiff)
    buildHeatmap(args.corrModule,outFileHeatmap)
    buildHistograms(outFileHistogram)



def buildClusterCaryotype(geneModuleFileName, outFileCluster, outFileCaryotype):
    """ From a list of gene-module association, this function creates two files needed by Circos, one keeping the gene-module association
        and adding simulated coordinates, the other one storing information about the modules."""
    global contenuCluster
    global contenuCaryotype
    
    geneModuleFile=open(geneModuleFileName,"r")

    prevModule = ""
    genePos = 0
    for currentLine in geneModuleFile.readlines()[1:]:
        
        ## Removal of extra characters
        currentLine = currentLine.strip()
        currentLine = currentLine.replace('"','')
        currentLine = currentLine.split(",")
        
        ## When reading next module, one needs to write information and set genePos back to 0
        if currentLine[1] != prevModule:
            if prevModule != "":
                contenuCaryotype += "chr - " + prevModule + " " + prevModule + " 0 " + str(genePos) + " " + prevModule + "\n"
            prevModule = currentLine[1]
            genePos = 0
        
        ## Get gene name and remove ENSEMBL Id if necessary    
        name = currentLine[0]
        pos = name.find("|")
        if (pos != -1):
            name = name[19:]        
        
        ## Store information about the gene
        contenuCluster += name+"\t"+currentLine[1]+"\t"+str(genePos)+"\t"+str(genePos + 1)+"\n"
        genePos += 1
    ## end for
    geneModuleFile.close()
    
    ## Store information about last module
    contenuCaryotype += "chr - " + currentLine[1]+ " " + currentLine[1] + " 0 " + str(genePos) + " " + currentLine[1]    
    
    ## Remove last carriage return and write to file
    contenuCluster = contenuCluster[:-1]
    with open(outFileCluster,"w") as fic:
        fic.write(contenuCluster)
        pass
    
    with open(outFileCaryotype,"w") as fic:
        fic.write(contenuCaryotype)
        pass





def buildGeneDiff(geneDiffFileName, outFileGeneDiff):
    """ From a list of differentially expressed genes, this function creates new file containing information about those genes."""
    global contenuCluster
    global contenuGeneDiff

    
    geneDiffFile=open(geneDiffFileName,"r")
    for currentGene in geneDiffFile.readlines():
        ## Removal of extra characters
        currentGene = currentGene.strip()
        currentGene = currentGene.replace('"','')
    
        ## Get gene name and remove ENSEMBL Id if necessary    
        pos = currentGene.find("|")
        if (pos != -1):
            currentGene=currentGene[19:]
            
        ## Get info for the gene diff    
        for i,geneInfo in enumerate(contenuCluster.split("\n")) :
            geneInfo = geneInfo.split("\t")
            if (geneInfo[0] == currentGene) :
                contenuGeneDiff += geneInfo[1] + " " + geneInfo[2] + " " + geneInfo[3] + " " + geneInfo[0] + "\n"
        
    geneDiffFile.close()    
    ## Remove last carriage return and write to file
    contenuGeneDiff = contenuGeneDiff[:-1]
    with open(outFileGeneDiff,"w") as fic:
        fic.write(contenuGeneDiff)
        pass    



def buildHeatmap(corrFileName, outFileHeatmap):
    """ This function creates a file containing information about the modules (including the correlation used for the heatmap). This will be used
        by Circos to display a heatmap track"""
    global contenuCaryotype

    contenuCor = ""
    correlation = open (corrFileName,"r")
    
    for curModule in correlation.readlines():
        ## Removal of extra characters
        curModule = curModule.strip()
        curModule = curModule.replace('\"','')
        curModule = curModule.replace('ME','')
        curModule = curModule.split(",")
        
        ## Get info about the modules
        for i,moduleInfo in enumerate(contenuCaryotype.split("\n")) :
            moduleInfo = moduleInfo.split(" ")
            if (moduleInfo[3] == curModule[0]) :
                    contenuCor += moduleInfo[3] + " " + moduleInfo[4] + " " + moduleInfo[5] + " " + curModule[1] + "\n"
    
    correlation.close()
                    
    ## Remove last carriage return and write to file
    contenuCor = contenuCor[:-1]                
    with open(outFileHeatmap,"w") as fic:
        fic.write(contenuCor)
        pass
                
                
def buildHistograms(outFileNameHistogram):
    """ This function computes the number and the proportion of differential genes in each module and generates two files
        needed by Circos to represent histograms """
    global contenuGeneDiff
    global contenuCaryotype
    
    nbTotGeneDiff = len(contenuGeneDiff.split("\n")) 
    
    histoCount = ""
    histoProp = ""
    for i,moduleInfo in enumerate(contenuCaryotype.split("\n")) :
        moduleInfo = moduleInfo.split(" ")
        nbGeneDiffInModule = 0
        for j,geneInfo in enumerate(contenuGeneDiff.split("\n")) :
            geneInfo = geneInfo.split(" ")
            if (geneInfo[0] == moduleInfo[3]):
                nbGeneDiffInModule += 1
        ##endfor
        ## Compute the ratio between DE genes in this module and total DE genes
        prop = float(nbGeneDiffInModule)/float(nbTotGeneDiff)
       
        histoCount += moduleInfo[3] + " " + moduleInfo[4] + " " + moduleInfo[5] + " " + str(nbGeneDiffInModule) + "\n"
        histoProp += moduleInfo[3] + " " + moduleInfo[4] + " " + moduleInfo[5] + " " + str(prop) + "\n"
    #endfor
    
    ## Remove last carriage return and write to file
    histoCount = histoCount[:-1]
    with open(outFileNameHistogram+"Count.txt","w") as fic:
        fic.write(histoCount)
        pass

    ## Remove last carriage return and write to file
    histoProp = histoProp[:-1]
    with open(outFileNameHistogram+"Prop.txt","w") as fic:
        fic.write(histoProp)
        pass

            
             
                
                
if __name__ == "__main__":
    main()
