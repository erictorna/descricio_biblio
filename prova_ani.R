library(dplyr)
library(RISmed)
library(easyPubMed)
library(rcrossref)
library(PubMedWordcloud)
library(wordcloud2)
library(parallel)
library(foreach)
library(doParallel)
library(tidytext)
library(tidyr)
library(ggplot2)
library(tm)
library(topicmodels)
library(viridis)
library(forcats)

# Primer definim la cerca tal com la posariem al PubMed
# search_topic <- '"Homosexuality, Female"[Mesh] OR lesbian* [ti] OR “women loving women”[tiab] OR “women who have sex with Women” [ti] OR WSW [ti] OR ((women [ti] OR female [ti]) AND (bicurious [tiab] OR heteroflexible [tiab] OR bisexual*[ti] OR bisexuality[MeSH Terms] OR homosexuality [MeSH Terms] OR homosexual*[ti] OR queer [ti] OR “same gender loving”[ti] OR “same sex attracted”[ti] OR “same sex couple”[ti] OR “same sex couples”[ti] OR “same sex relations”[ti] OR “sexual and gender minorities”[ti] OR “sexual and gender minority”[ti] OR “sexual identity”[ti] OR “sexual minorities”[ti] OR “sexual minority”[ti] OR “sexual orientation”[ti])) NOT (gay  [ti] OR MSM [ti] OR “men who have sex with men” [ti]) '
search_topic <- ' (("homosexuality, female"[MeSH Terms] OR "lesbian*"[Title] OR "women loving women"[Title/Abstract] OR "women who have sex with Women"[Title] OR "WSW"[Title] OR (("women"[Title] OR "female"[Title]) AND ("bicurious"[Title/Abstract] OR "heteroflexible"[Title/Abstract] OR "bisexual*"[Title] OR "bisexuality"[MeSH Terms] OR "homosexuality"[MeSH Terms] OR "homosexual*"[Title] OR "queer"[Title] OR "same gender loving"[Title] OR "same sex attracted"[Title] OR "same sex couple"[Title] OR "same sex couples"[Title] OR "same sex relations"[Title] OR "sexual and gender minorities"[Title] OR "sexual and gender minority"[Title] OR "sexual identity"[Title] OR "sexual minorities"[Title] OR "sexual minority"[Title] OR "sexual orientation"[Title]))) NOT ("gay"[Title] OR "MSM"[Title] OR "men who have sex with"[Title] OR "transgender women"[Title] OR "transgender men"[Title] OR "trans men"[Title] OR "trans women"[All Fields] OR "lgbt*"[Title] OR "homosexuality, male"[MeSH Terms])) AND ((humans[Filter]) AND (2004:3000/12/12[pdat]) AND (english[Filter])) '
my_query <- get_pubmed_ids(search_topic)
my_abstracts_xml <- fetch_pubmed_data(my_query,retmax = as.numeric(my_query$Count))
all_xml <- articles_to_list(my_abstracts_xml,simplify = F)
search_query <- EUtilsSummary(search_topic, retmax=as.numeric(my_query$Count))
records <- EUtilsGet(search_query, type = "efetch", db = "pubmed")
# ---------------------

# función para extraer el número de citaciones por DOI en rcrossref sin errores
Ncitation <- function(doi){
  if(nchar(doi)==0){
    cit <- data.frame(doi=doi,count=NA)
  }else{
    cit <- cr_citation_count(doi,async = T)
  }
  names(cit) <- c("doi","Num_Citation")
  return(cit)
}

pubmed_data <- data.frame(pmid=PMID(records),
                          Language = Language(records),
                          country=Country(records),stringsAsFactors = F)

# Proceso de paralelización para construir la base de datos final
cl <- makeCluster(8) 
registerDoParallel(cl)

fullDF <- tryCatch(
  {foreach(x=all_xml, 
           .packages = 'easyPubMed',
           .combine = rbind) %dopar% article_to_df(
             pubmedArticle = x,
             autofill = T, 
             max_chars = -1, #Extrae el abstract completo
             getKeywords = T,
             getAuthors = F)},
  error = function(e) {NULL},
  finally = {stopCluster(cl)}) 

# Extraer información del número de citaciones del DOI
cit <- do.call("rbind",lapply(unique(fullDF$doi),Ncitation))


tyPub <- data.frame(pmid=records@PMID,
                    PubType=sapply(records@PublicationType,
                                   function(x) paste(x,collapse = " ")),
                    stringsAsFactors = F)


FinalDF <- left_join(x = fullDF,y = pubmed_data,"pmid")
FinalDF <- left_join(x = FinalDF,y = cit,"doi")
FinalDF <- left_join(x = FinalDF,y = tyPub,"pmid")

# Terminos Mesh
Mesh <- records@Mesh
names(Mesh) <- records@PMID

rm(my_abstracts_xml,all_xml)
rm(cit,fullDF,pubmed_data)

FinalDF %>% mutate(abstract=substr(abstract,1,40),
                   title=substr(title,1,20),
                   keywords=substr(keywords,1,20)) %>%
  distinct(pmid,.keep_all = T) %>% DT::datatable()
#----------------------------
# Núvol de paraules amb Mesh més utilitzats
library(data.table)
BDMesh <- rbindlist(Mesh, idcol = names(Mesh)) %>%
  `colnames<-`(c("pmid", "Heading", "Type"))

Qualifier <- BDMesh %>% filter(Type=="Qualifier") %>%
  dplyr::select(Heading) %>% unlist %>% as.character  %>% cleanAbstracts
wordcloud2(Qualifier,size=1.5,color='random-dark',gridSize = 10)
#-----------------------------------
Meshtable <- BDMesh[, .(Heading = toString(Heading)), by = pmid]
FinalDF = left_join(FinalDF, Meshtable)
FinalDF = FinalDF %>% rename(Mesh = Heading)
# Creem una taula que agrupi només els articles que incloguin certa paraula, o grup d'elles, en l'abstract, títol, Mesh, etc
FinalDF[grepl("Male", FinalDF$Mesh),] %>%
  distinct(pmid,.keep_all = T) %>%
  dplyr::select(title, journal, doi, Mesh) %>%
  DT::datatable(options = list(scrollX = TRUE,
                               initComplete = DT::JS(
                                 "function(settings, json) {",
                                 "$(this.api().table().container()).css({'font-size': '70%'});",
                                 "}")))

nomes_mesh = FinalDF %>% dplyr::select(title, doi, Mesh)
rm(cl, Meshtable, my_query, records, search_query, tyPub, Mesh)
# Obtenim un taula amb tots els mesh existents

library(xml2)
library(XML)
library(tibble)
library(stringr)
# https://www.nlm.nih.gov/mesh/intro_record_types.html
meshes = read_xml('~/idiap/projects/descrip_pubmed/desc2024.xml')
meshes_xml <- xmlParse(meshes)
data = xmlToDataFrame(nodes = getNodeSet(meshes_xml, "//DescriptorRecord"))
data = data %>% dplyr::select(DescriptorName, TreeNumberList)
setDT(data)
C = data[like(TreeNumberList, "C")]
E = data[like(TreeNumberList, "E")]
N = data[like(TreeNumberList, "N")]
setDT(nomes_mesh)

nomes_mesh[, Mesh := strsplit(Mesh, ",")]

nomes_mesh_expanded <- nomes_mesh[, .(doi, Mesh = unlist(Mesh)), by = doi]
nomes_mesh_expanded[, Mesh := trimws(Mesh)]

# setkey(nomes_mesh_expanded, Mesh)
# setkey(diseases, DescriptorName)

result_C = nomes_mesh_expanded[Mesh %in% C$DescriptorName]
result_E = nomes_mesh_expanded[Mesh %in% E$DescriptorName]
result_N = nomes_mesh_expanded[Mesh %in% N$DescriptorName]


fwrite(result_C, file = '~/idiap/projects/descrip_pubmed/diseases.csv')
fwrite(result_E, file = '~/idiap/projects/descrip_pubmed/Analytical_diagnostic_techiques.csv')
fwrite(result_N, file = '~/idiap/projects/descrip_pubmed/HealthCare.csv')
