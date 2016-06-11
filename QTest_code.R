##~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#  R code for QTest            
#                                
#  by Jaehoon Lee              
#                              
#   jhlee1213@gmail.com        
##~~~~~~~~~~~~~~~~~~~~~~~~~~~####

library(plyr)
library(CCP)
library(homals)
library(Matrix)
library(SKAT)
library(CompQuadForm)
library(parallel)

## Empirical null distribution for Q3 test
null.dist.Q3<-read.table("n.log10.minp.1e09.txt",header=T)

## Get a given STT
get.a<-function(L,STT){
	aa<-diff<-seq(0,1,length=200)
	for(i in 1:length(aa)){
	diff[i]<-abs(get.stt(L,aa[i],STT)-STT)
	}
	return(aa[which.min(diff)])
}
get.stt<-function(L,a,STT){
	1-pgamma(L*qgamma(1-STT,a,1),L*a,1)
}
## get minor allele frequency
get.maf <-function(vec){
	(length(which(vec==1))+2*length(which(vec==2)))/(2*length(vec))
}

###############################################################################################
## QTest.one(phe.cova,geno,yname, r2)  #QTest for single gene                                
##                                                                                           
## y: trait, newgeno: genotype matrix                                           
## cut.r2: r^2 cutoff, n.perm: # of pumutations                                              
## n.perm: number of permutation for GM method                                              
## a: shape parameter for GM method                                                          
###############################################################################################

QTest.one<-function(y,covadat=NULL,newgeno,STT,weight=FALSE){
	if(length(covadat)!=0){resid<-try(resid(glm(y~.,data=data.frame(covadat),na.action=na.exclude)),TRUE)}
	if(length(covadat)==0){resid<-try(resid(glm(y~1,na.action=na.exclude)),TRUE)}
	S<-try(as.matrix(newgeno),TRUE)
	fit<-try(glm(resid~.,data=data.frame(S)),TRUE)
	na.S<-try(which(is.na(coef(fit)[-1])==TRUE),TRUE)
	if(length(na.S)>0){
		S<-try(as.matrix(S[,-na.S]),TRUE)
		fit<-try(glm(resid~.,data=data.frame(S)),TRUE)
	}
	coef<-try(coef(summary(fit))[-1,1:2],TRUE)
	
	
	if(mode(fit)=="character"){length(coef)<-0}
	if(length(coef)!=0){
		if(length(coef)!=2){beta1<-coef[,1];se1<-coef[,2]}
		if(length(coef)==2){beta1<-coef[1];se1<-coef[2]}
		SS<-cbind(1,S); n<-length(y)
		vv<-vcov(fit)[-1,-1];alpha<-(1/(se1^2))
		#vv<-solve(t(SS)%*%SS)[-1,-1]*var(y)*(n-1)/n;alpha<-(1/(se1^2))
	
		##QTest1##
		if(weight==FALSE){
			var.pool0<-t(alpha)%*%vv%*%alpha
			beta.pool0<-t(alpha)%*%beta1
			z.score0<-beta.pool0/sqrt(var.pool0)
			Q1<-z.score0^2
			p.Q1<-pchisq(Q1,df=1,lower.tail=FALSE)
		}
		if(weight==TRUE){
			maf.S<-apply(S,2,function(v)mean(na.omit(v))/2)
			#w.S0<-1/sqrt(maf.S*(1-maf.S));
			w.S0<-qbeta(maf.S,1,25,lower.tail=F)
			WS<-diag(w.S0);if(length(beta1)==1){WS<-w.S0}
			var.pool<-t(alpha)%*%WS%*%vv%*%WS%*%alpha
			beta.pool<-t(alpha)%*%WS%*%beta1
			z.score<-beta.pool/sqrt(var.pool)
			Q1<-z.score^2
			p.Q1<-pchisq(Q1,df=1,lower.tail=FALSE)
		}
		
		## QTest2 ##
		Q2.eigen<-eigen(vv);U2<-Q2.eigen$vectors; l2<-Q2.eigen$values;na.l2<-which(l2/mean(l2)<0.01)
		if(length(na.l2)>0){l2<-l2[-na.l2];U2<-U2[,-na.l2]}
		p2<-pchisq((t(U2)%*%beta1)^2/l2,df=1,lower.tail=FALSE)
		a<-get.a(length(p2),STT)
		q2<-2*(qgamma(p2,a,1,lower.tail=FALSE))
		Q2<-sum(q2)
		p.Q2<-pchisq(Q2,df=2*a*length(q2),lower.tail=FALSE)


		##QTest3 ##

		if(weight==FALSE){
			c(t(alpha)%*%vv)->cov.beta
			b.star<-beta1-beta.pool0*cov.beta/var.pool0[1]
			vv.star<-vv-(cov.beta%*%t(cov.beta))/var.pool0[1]
		}
		if(weight==TRUE){
			w.vv<-WS%*%vv%*%WS
			c(t(alpha)%*%w.vv)->cov.beta
			b.star<-WS%*%beta1-beta.pool*cov.beta/var.pool[1]
			vv.star<-w.vv-(cov.beta%*%t(cov.beta))/var.pool[1]
		}

		if(length(beta1)!=1){
			Q3.eigen<-eigen(vv.star); U3<-Q3.eigen$vectors; l3<-Q3.eigen$values;na.l3<-which(l3/mean(l3)<0.001)
			if(length(na.l3)>0){l3<-l3[-na.l3];U3<-U3[,-na.l3]}
			L3<-try(diag(l3),TRUE);if(length(l3)==1){L3<-l3}
			q2.proj<-(t(U3)%*%b.star)^2/l3; 
			p2.1<-pchisq(q2.proj,df=1,lower.tail=FALSE)
			a<-get.a(length(p2.1),STT)
			Q2.proj<-sum(2*qgamma(p2.1,a,1,lower.tail=FALSE))
			p.Q2.proj<-pchisq(Q2.proj,df=2*a*length(l3),lower.tail=FALSE)
			Q2.1<-qchisq(p.Q2.proj,df=1,lower.tail=FALSE)
		
			pi0<-seq(0,1,by=0.1);p.Q3.can<-Q3<-rep(1,11)
			p.Q3.can[1]<-p.Q2.proj; Q3[1] <- Q2.1
			p.Q3.can[11]<-p.Q1; Q3[11]<-Q1
			for(h in 2:10){
				Q3[h]<-pi0[h]*Q1+(1-pi0[h])*Q2.1
				p.Q3.can[h]<-davies(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))$Qq
				if(p.Q3.can[h]<=0|p.Q3.can[h]>1){p.Q3.can[h]<-imhof(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))$Qq}
				if(p.Q3.can[h]<=0|p.Q3.can[h]>1){p.Q3.can[h]<-liu(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))[1]}
			}
	
	
			Q3final<-Q3[which.min(p.Q3.can)]
			p.Q3<-(sum(null.dist.Q3[null.dist.Q3[,1]>-log10(min(p.Q3.can)),2])+1)/(sum(null.dist.Q3[,2])+1)
	
		}

		if(length(beta1)==1){p.Q3<-p.Q1; Q3final<-Q1}
		rslt<-list( data.frame(Q1, Q2, Q3=Q3final), data.frame(p.Q1,p.Q2,p.Q3))
		names(rslt)<-c("Qstatistic", "p.value")	
	}

	if(length(coef)==0){rslt<-NA}
	options (warn=-1)
	return(rslt)	
}



SKAT.one<-function(y,covadat=NULL,geno,weights.beta=c(1,25)){
	if(length(covadat)>0){obj.jh<-SKAT_Null_Model(y~.,data=data.frame(covadat),out_type="C")}
	if(length(covadat)==0){obj.jh<-SKAT_Null_Model(y~1,out_type="C")}
	p.skat<-try(SKAT(as.matrix(geno), obj.jh,weights.beta=weights.beta)$p.value,TRUE)
	p.skat.o<-try(SKAT(as.matrix(geno), obj.jh,weights.beta=weights.beta,method="optimal")$p.value,TRUE)
	return(c(p.skat,p.skat.o))
}


###################################
##    Gene-SNP annotation file                   
##    column names:  "CHR" "gene","snp" 
###################################
#anno<-read.table("...",header=T)

###################################################################
##    Pheno & covariate file                                
##   "ID": 1st column for sample ID 
##    column 2nd, 3rd, ... : phenotype or covariate 
###################################################################
#phe.cova<-read.table("...",header=T)

########################################
##    Genotype file for all SNPs      ##  
##    "ID": 1st column for sampleID
##    column 2nd, 3rd, ...: SNP genotypes
########################################
#allgeno<-read.table("...",header=T)

#######################################
## Pathway-Gene annotation file
## column names: "setName", "gene1", "gene2", ...
#######################################

###############################################################################################
## QTest for gene set analysis     (GSQ.QTest and GSB.QTest)                      
##                                                                                           
## yname: trait name in phe.cova file                                          
## covaname: covariate name in phe.cova file
## cut.r2: r^2 cutoff, n.perm: # of pumutations                                              
## genePerm: whether permutation is used for gene-level
## setPerm: whether permutation is used for set-level
## nPerm: number of permutation for GM method                                              
## a: shape parameter for GM method                                                          
## method: "BQ": GSB.QTest, "QQ": GSQ.QTest
## omit.onesizegene: whether omit signle size genes in analysis
## outname: name for output file
## mc.cores: number of multi-cores for analysis
###############################################################################################

GS.QTest<-function(phe.cova,yname,covaname=NULL,allgeno,set.anno,gene.anno,maxSetSize=200,STT=0.2,maf.cut=0.05,min.mac=5,weight=FALSE,method=c("QQ","BQ"),weight.type="inv.var",omit.onesizegene=FALSE, genePerm=FALSE, setPerm=FALSE, mc.cores=10, nPerm, outname){
	cat("spliting genotype data for each gene..","\n")
	allgeno<-allgeno[match(phe.cova$ID,allgeno$ID),]
	allgeno<-allgeno[,-1]
	where.y<-which(colnames(phe.cova)==yname)
	where.cova<-which(colnames(phe.cova)%in%covaname==TRUE)
	y<-phe.cova[,where.y]
	covadat<-phe.cova[,where.cova]
	na.y<-which(is.na(y)==TRUE)	
	if(length(na.y)>0){
		y<-y[-na.y]
		if(length(covaname)!=0){covadat<-covadat[-na.y,]}
		allgeno<-allgeno[-na.y,]
	}
	allgeno<-allgeno[,apply(allgeno,2,function(v)length(which(v>0)))>=min.mac]
	if(length(covaname)==0){covadat<-NULL}
	frq<-apply(allgeno,2,get.maf)
	allgeno<-allgeno[,frq<maf.cut]
	gene.anno<-gene.anno[gene.anno$snp%in%colnames(allgeno)==TRUE,]
	if(omit.onesizegene==TRUE){
		gene.anno$gene<-as.character(gene.anno$gene);gene.anno$snp<-as.character(gene.anno$snp);gene.anno$CHR<-as.character(gene.anno$CHR)
		na.gene<-names(table(gene.anno$gene)[table(gene.anno$gene)==1])
		gene.anno<-gene.anno[gene.anno$gene%in%na.gene==FALSE,]
	}
	geno.set1<-tapply(gene.anno$snp,gene.anno$gene,function(v)allgeno[,colnames(allgeno)%in%v==TRUE])	
	genename<-names(geno.set1)

	rslt.set<-list()
	apply(set.anno[,-1],1,function(v)na.omit(as.character(as.matrix(v))))->set.list
	if(class(set.list)!="list"){apply(set.anno[,-1],1,function(v)na.omit(as.character(as.matrix(v))))->set.list}
	names(set.list)<-set.anno[,1]
	set.list1<-lapply(set.list,function(v)v[v%in%genename==TRUE])
	geneleng<-unlist(lapply(set.list1,length))
	geneleng1<-geneleng[match(set.anno[,1],names(geneleng))];geneleng1[is.na(geneleng1)==TRUE]<-0
	idx.accept<-which(geneleng1<maxSetSize&geneleng1>1)

	for(i in idx.accept){
		if(method=="QQ"){
			conv.idx<-which(names(set.list1)%in%names(geneleng1)[i]==TRUE)
			setgeno<-apply(as.matrix(set.list1[[conv.idx]]),1,function(v)list(as.matrix(geno.set1[names(geno.set1)%in%v==TRUE][[1]])))
			setgeno<-lapply(setgeno,function(v)v[[1]])
			names(setgeno)<-data.frame(set.list1[[conv.idx]])[,1]
			setleng0<-unlist(lapply(setgeno,ncol))
			idx.l1<-which(setleng0==1)
			if(length(idx.l1)>0){for(k in 1:length(idx.l1)){colnames(setgeno[[idx.l1[k]]])<-gene.anno$snp[gene.anno$gene==names(setgeno)[idx.l1[k]]]}}

			naIdx<-list()
			hmatList<-list()
			hmatSample<-list()
			hmatSample1<-list()

			for(j in 1:length(setgeno)){
				setgeno[[j]]<-as.matrix(setgeno[[j]])
				cmat<-as.matrix(setgeno[[j]])		
				n.col<-ncol(setgeno[[j]]); n.row<-nrow(setgeno[[j]])
				cmatSample<-1:n.row
				naSampleIdx<-which(apply(cmat,1,function(v)length(which(is.na(v)==TRUE)))==ncol(cmat))
				hmatSample[[j]]<-cmatSample
				if(length(naSampleIdx)>0){
					if(length(naSampleIdx)==n.row){
						naIdx[[j]]<-TRUE
					} else {
						naIdx[[j]]<-FALSE
						setgeno[[j]]<-as.matrix(setgeno[[j]][-naSampleIdx,])
						hmatSample[[j]]<-cmatSample[-naSampleIdx]
					}
				} else {naIdx[[j]]<-FALSE}
			}
			names(hmatSample)<-names(setgeno)
			if(length(which(naIdx==TRUE))>0){
				setgeno<-setgeno[-which(naIdx==TRUE)]
				hmatSample<-hmatSample[-which(naIdx==TRUE)]
			}

			for(j in 1:length(setgeno)){
				cmat<-as.matrix(setgeno[[j]])
				idx.leng0<-which(apply(cmat,2,sum)==0);if(length(idx.leng0)>0){cmat<-cmat[,-idx.leng0]}
				if(ncol(cmat)==1){
					options(warn=-1)
					capture.output(hmat <- homals(data.frame(cmat), eps=1e-05, ndim=1)$obj, file="NUL")
				} else {
					options(warn=-1)
					capture.output(hmat <- homals(data.frame(cmat), eps=1e-05)$obj, file="NUL")
				}
				hmatList[[j]]<-hmat
			}
			names(hmatList)<-names(setgeno)
			comb<-combn(length(setgeno),2)
			ncomb<-ncol(comb)
			p.comb<-rep(NA,ncomb)

			for(j in 1:ncomb){
				hmat1<-hmatList[[comb[,j][1]]]; 
				hmat1Sample<-hmatSample[[comb[,j][1]]]
				hmat2<-hmatList[[comb[,j][2]]]; 
				hmat2Sample<-hmatSample[[comb[,j][2]]]
			
				shareSample<-intersect(hmat1Sample, hmat2Sample)
				hmat1<-as.matrix(hmat1[match(shareSample, hmat1Sample),])
				hmat2<-as.matrix(hmat2[match(shareSample, hmat2Sample),])

				can.cor<-try(cancor(hmat1, hmat2)$cor, TRUE)
				if(mode(can.cor)=="character"){
					print("NA cancor");
				} else {
					if(length(can.cor)<min(ncol(hmat1),ncol(hmat2))){can.cor<-c(can.cor,rep(0,min(ncol(hmat1),ncol(hmat2))-length(can.cor)))}
					capture.output(asym<-p.asym(rho=can.cor, nrow(hmat1), ncol(hmat1),ncol(hmat2),tstat="Wilks"), file="NUL")
					p.comb[j]<-pf(asym$approx,asym$df1,asym$df2,lower.tail=F)[1]
				}

			}

			comb1<-apply(as.matrix(comb[,p.comb<0.05/ncomb]),2,function(v)list(na.omit(v)))
			comb1<-lapply(comb1,unlist)
			if(length(which(lapply(comb1,length)==0))>0){comb1<-comb1[-which(lapply(comb1,length)==0)]}
	
			comb.tab<-table(unlist(comb1))
			comb.replicate<-as.numeric(names(comb.tab)[comb.tab>1])
			if(length(comb.replicate)>0){
				for(j in 1:length(comb.replicate)){
					comb1.idx=which(lapply(comb1,function(v)length(which(v%in%comb.replicate[j]==TRUE)))>0)
					temp.comb1<-comb1[comb1.idx]
					comb1<-comb1[-comb1.idx]
					comb1[[length(comb1)+1]]<-unique(unlist(temp.comb1))
				}
			}

			options(warn=-1)
			unique.comb<-na.omit(unique(unlist(comb1)))
			if(length(unique.comb)>0){	
				setgeno1<-setgeno[-unique.comb]
				hmatSample1<-hmatSample[-unique.comb]
				for(j in 1:length(comb1)){
					intersec<-1:length(y)
					for(k in 1:length(comb1[[j]])){
						intersec<-intersect(intersec,hmatSample[[comb1[[j]][k]]])
					}
					tempgeno<-list()
					for(k in 1:length(comb1[[j]])){
						tempgeno[[k]]<-as.matrix(setgeno[[comb1[[j]][k]]])[match(intersec,hmatSample[[comb1[[j]][k]]]),]
								}
					
					setgeno1[[length(setgeno)-length(unique.comb)+j]]<-do.call(cbind,tempgeno)
					hmatSample1[[length(setgeno)-length(unique.comb)+j]]<-intersec
					names(setgeno1)[length(setgeno)-length(unique.comb)+j]<-paste(names(setgeno)[comb1[[j]]],collapse="_")
				}	
			} else {
				hmatSample1<-hmatSample
				setgeno1<-setgeno
			}

			if(genePerm==TRUE){
				rsltP0<-mclapply(1:length(setgeno1),function(k)QTestPerm.one(y=y,newgeno=as.matrix(setgeno1[[k]]),nPerm, STT=STT,weight=weight,mc.cores=mc.cores), mc.cores=1)
				rsltP<-do.call(rbind,rsltP0)
				write.table(rsltP,paste(names(idx.accept)[i],"_permuteResult.txt",sep=""),quote=F,row.names=F)
			}else{
				rslt<-t(mapply(function(v1,v2)QTest.one(y=y[v2],newgeno=as.matrix(v1),STT=STT,weight=weight)$p.value, setgeno1,hmatSample1))
				rslt.skat<-t(mapply(function(v1,v2)SKAT.one(y=y[v2],geno=as.matrix(v1)), setgeno1,hmatSample1))
			}
			
			a<-get.a(length(setgeno1),STT)
			q1<-2*(qgamma(as.numeric(rslt[,1]),a,1,lower.tail=FALSE))
			q2<-2*(qgamma(as.numeric(rslt[,2]),a,1,lower.tail=FALSE))
			q3<-2*(qgamma(as.numeric(rslt[,3]),a,1,lower.tail=FALSE))
			q.skat<-2*(qgamma(as.numeric(rslt.skat[,1]),a,1,lower.tail=FALSE))
			q.skato<-2*(qgamma(as.numeric(rslt.skat[,2]),a,1,lower.tail=FALSE))
	
			Q1<-sum(q1);Q2<-sum(q2);Q3<-sum(q3)
			Q.skat<-sum(q.skat); Q.skato<-sum(q.skato)
			p.Q1<-pchisq(Q1,df=2*a*length(q1),lower.tail=FALSE)
			p.Q2<-pchisq(Q2,df=2*a*length(q1),lower.tail=FALSE)
			p.Q3<-pchisq(Q3,df=2*a*length(q1),lower.tail=FALSE)
			p.skat<-pchisq(Q.skat,df=2*a*length(q.skat),lower.tail=FALSE)
			p.skato<-pchisq(Q.skato,df=2*a*length(q.skato),lower.tail=FALSE)
			rslt.set[[i]]<-data.frame(set=names(set.list1[conv.idx]), n.gene=length(setgeno1), n.SNP=sum(unlist(lapply(setgeno1,ncol))),p.Q1,p.Q2,p.Q3, p.skat,p.skato)
			if(setPerm==TRUE){
				rsltPP <- mclapply(1:nPerm, function(k){
				yP<-y[sample(1:length(y),length(y),replace=F)];
				rsltP<-do.call(rbind,lapply(setgeno1,function(v)QTest.one(y=yP,newgeno=as.matrix(v),STT=STT,weight=weight)$p.value))
				a<-get.a(length(setgeno1),STT)
				q1P<-2*(qgamma(rsltP[,1],a,1,lower.tail=FALSE))
				q2P<-2*(qgamma(rsltP[,2],a,1,lower.tail=FALSE))
				q3P<-2*(qgamma(rsltP[,3],a,1,lower.tail=FALSE))	
				Q1P<-sum(q1P);Q2P<-sum(q2P);Q3P<-sum(q3P)
				c(Q1P, Q2P, Q3P)
			}, mc.cores=mc.cores)
			rsltPP1<-do.call(rbind,rsltPP)
			p.Q1P<-length(which(rsltPP1[,1]>Q1))/nrow(rsltPP1)
			p.Q2P<-length(which(rsltPP1[,2]>Q2))/nrow(rsltPP1)
			p.Q3P<-length(which(rsltPP1[,3]>Q3))/nrow(rsltPP1)
			rslt.set[[i]]<-data.frame(set=names(set.list1[conv.idx]), n.gene=length(setgeno1), n.SNP=ncol(do.call(cbind,setgeno1)),p.Q1,p.Q2,p.Q3, p.Q1P, p.Q2P, p.Q3P)
		}
		
	}
	
	if(method=="BQ"){
		conv.idx<-which(names(set.list1)%in%names(geneleng1)[i]==TRUE)	
		setgeno<-apply(as.matrix(set.list1[[conv.idx]]),1,function(v)list(as.matrix(geno.set1[names(geno.set1)%in%v==TRUE][[1]])))
		setgeno<-lapply(setgeno,function(v)v[[1]])
		names(setgeno)<-data.frame(set.list1[[conv.idx]])[,1]
		setleng0<-unlist(lapply(setgeno, ncol))
		idx.l1<-which(setleng0==1)
		if(length(idx.l1)>0){for(k in 1:length(idx.l1)){colnames(setgeno[[idx.l1[k]]])<-gene.anno$snp[gene.anno$gene==names(setgeno)[idx.l1[k]]]}}
		macS<-lapply(setgeno,function(v0)apply(as.matrix(v0),2,function(v)length(which(na.omit(v)>0))))
		naS<-lapply(macS, function(v)which(v<min.mac))
		set.snpname<-lapply(setgeno,colnames)
		setgeno2<-setgeno1<-mapply(function(v1,v2){if(length(v2)>0)return(v1[,-v2]);if(length(v2)==0)return(v1)},setgeno,naS,SIMPLIFY=F)
		set.snpname1<-mapply(function(v1,v2){if(length(v2)>0)return(v1[-v2]);if(length(v2)==0)return(v1)},set.snpname,naS)
		if(mode(setgeno2)!="list"){	setgeno2<-setgeno1<-mapply(function(v1,v2){if(length(v2)>0)return(list(v1[,-v2]));if(length(v2)==0)return(v1)},setgeno,naS)}
		if(mode(set.snpname1)!="list"){	set.snpname1<-mapply(function(v1,v2){if(length(v2)>0)return(list(v1[-v2]));if(length(v2)==0)return(v1)},set.snpname,naS)}
		setleng1<-lapply(setgeno2,function(v)ncol(as.matrix(v)))
		if(length(which(setleng1==0))>0){
			setgeno2<-setgeno1[-which(setleng1==0)]
			set.snpname1<-set.snpname1[-which(setleng1==0)]
			setleng1<-lapply(setgeno2,function(v)ncol(as.matrix(v)))
		}
		if(length(setgeno2)>0){	
			if(weight==TRUE){
				mafS1<-lapply(setgeno2,function(v0)apply(as.matrix(v0),2,function(v)mean(na.omit(v))/2))
				if(weight.type=="inv.var"){wS<-lapply(mafS1,function(v)1/sqrt(v*(1-v)))}
				if(weight.type=="beta"){wS<-lapply(mafS1,function(v)qbeta(v,1,25,lower.tail=F))}
				wS<-lapply(wS,function(v)v/sum(v))
				colS<-mapply(function(v1,v2)apply(as.matrix(v1),1,function(v0)sum(v2*v0)),setgeno2,wS)
			}
			if(weight==FALSE){colS<-do.call(cbind,lapply(setgeno2,function(v1)apply(as.matrix(v1),1,function(v0)sum(v0))))}

			if(setPerm==FALSE){
				rslt<-try(QTest.one(y,covadat=NULL,newgeno=colS,STT)$p.value,TRUE)
				rslt.skat<-try(SKAT.one(y,covadat=NULL,geno=colS),TRUE)
				if(mode(rslt)=="character"){rslt<-NA};if(mode(rslt.skat)=="character"){rslt.skat<-NA}
				if(is.na(rslt)==TRUE){rslt<-matrix(c(NA,NA,NA),1,3);colnames(rslt)<-c("p.Q1","p.Q2","p.Q3")}
				if(is.na(rslt.skat)==TRUE){rslt<-matrix(c(NA,NA),1,2);colnames(rslt)<-c("p.skat","p.skato")}
			} else {
				rsltPP <- do.call(rbind,mclapply(1:nPerm, function(k){
					yP<-y[sample(1:length(y),length(y),replace=F)]
					QTest.one(y=yP,covadat=NULL, newgeno=colS,STT=STT,weight=weight)$p.value
				},mc.cores=mc.cores))
				rslt0<-try(as.numeric(QTest.one(y,covadat=NULL,newgeno=colS,STT=STT,weight=weight)$p.value),TRUE)
				rslt <- data.frame(p.Q1=(length(which(rsltPP[,1]<rslt0[1]))+1)/(nPerm+1),
				           p.Q2=(length(which(rsltPP[,2]<rslt0[2]))+1)/(nPerm+1),
				           p.Q3= (length(which(rsltPP[,3]<rslt0[3]))+1)/(nPerm+1))
				}
				rslt.set[[i]]<-data.frame(set=names(set.list1[conv.idx]),n.gene=length(setgeno2),n.SNP=sum(unlist(setleng1)),rslt, p.skat=rslt.skat[1],p.skato=rslt.skat[2])
			}	
		}

		cat("Test for set",i,"\n");flush.console()
	}

	rslt.set1<-do.call(rbind,rslt.set)
	write.table(rslt.set1,paste("GS.QTest_",outname,".rslt",sep=""),quote=F,row.names=F)
	return(rslt.set1)
}

#Funtion  for identify signifcant factors within gene-set    

QTest.within.geneset<-function(phe.cova,yname,covaname=NULL,allgeno,set.anno,gene.anno,maxSetSize,STT,maf.cut,min.mac,weight=FALSE, weight.type="inv.var", omit.onesizegene=TRUE, genePerm=FALSE, nPerm, mc.cores, outname){
	cat("spliting genotype data for each gene..","\n")
	allgeno<-allgeno[match(phe.cova$ID,allgeno$ID),]
	allgeno<-allgeno[,-1]
	where.y<-which(colnames(phe.cova)==yname)
	where.cova<-which(colnames(phe.cova)%in%covaname==TRUE)
	y<-phe.cova[,where.y]
	covadat<-phe.cova[,where.cova]
	na.y<-which(is.na(y)==TRUE)	
	if(length(na.y)>0){
		y<-y[-na.y]
		if(length(covaname)!=0){covadat<-covadat[-na.y,]}
		allgeno<-allgeno[-na.y,]
	}
	allgeno<-allgeno[,apply(allgeno,2,function(v)length(which(v>0)))>=min.mac]
	if(length(covaname)==0){covadat<-NULL}
	frq<-apply(allgeno,2,get.maf)
	allgeno<-allgeno[,frq<maf.cut]
	gene.anno<-gene.anno[gene.anno$snp%in%colnames(allgeno)==TRUE,]
	gene.anno$gene<-as.character(gene.anno$gene);gene.anno$snp<-as.character(gene.anno$snp);gene.anno$CHR<-as.character(gene.anno$CHR)
	if(omit.onesizegene==TRUE){
		na.gene<-names(table(gene.anno$gene)[table(gene.anno$gene)==1])
		gene.anno<-gene.anno[gene.anno$gene%in%na.gene==FALSE,]
	}
	geno.set1<-tapply(gene.anno$snp,gene.anno$gene,function(v)allgeno[,colnames(allgeno)%in%v==TRUE])	
	genename<-names(geno.set1)

	rslt.set.q<-rslt.set.b<-rslt.gene<-list()
	cat("Preprocessing...","\n");flush.console();
	apply(set.anno[,-1],1,function(v)na.omit(as.character(unlist(v))))->set.list
	if(nrow(set.anno)==1){set.list<-list(set.list)}
	names(set.list)<-set.anno[,1]
	set.list1<-lapply(set.list,function(v)v[v%in%genename==TRUE])
	geneleng<-unlist(lapply(set.list1,length))
	geneleng1<-geneleng[match(set.anno[,1],names(geneleng))];geneleng1[is.na(geneleng1)==TRUE]<-0
	idx.accept<-which(geneleng1<maxSetSize&geneleng1>1)

	for(i in idx.accept){

		## method=="QQ"
		cat("Test for set",i,"\n");flush.console()
		conv.idx<-which(names(set.list1)%in%names(geneleng1)[i]==TRUE)
		setgeno<-apply(as.matrix(set.list1[[conv.idx]]),1,function(v)list(as.matrix(geno.set1[names(geno.set1)%in%v==TRUE][[1]])))
		setgeno<-lapply(setgeno,function(v)v[[1]])
		names(setgeno)<-data.frame(set.list1[[conv.idx]])[,1]

		naIdx<-list()
		hmatList<-list()
		hmatSample<-list()
		hmatSample1<-list()

		for(j in 1:length(setgeno)){
			setgeno[[j]]<-as.matrix(setgeno[[j]])
			cmat<-as.matrix(setgeno[[j]])		
			n.col<-ncol(setgeno[[j]]); n.row<-nrow(setgeno[[j]])
			cmatSample<-1:n.row
			naSampleIdx<-which(apply(cmat,1,function(v)length(which(is.na(v)==TRUE)))==ncol(cmat))
			hmatSample[[j]]<-cmatSample
			if(length(naSampleIdx)>0){
				if(length(naSampleIdx)==n.row){
					naIdx[[j]]<-TRUE
				} else {
					naIdx[[j]]<-FALSE
					setgeno[[j]]<-as.matrix(setgeno[[j]][-naSampleIdx,])
					hmatSample[[j]]<-cmatSample[-naSampleIdx]
				}
			} else {naIdx[[j]]<-FALSE}
		}
		setleng0<-unlist(lapply(setgeno,ncol))
		idx.l1<-which(setleng0==1)
		if(length(idx.l1)>0){for(k in 1:length(idx.l1)){colnames(setgeno[[idx.l1[k]]])<-gene.anno$snp[gene.anno$gene==names(setgeno)[idx.l1[k]]]}}

		names(hmatSample)<-names(setgeno)

		if(length(which(naIdx==TRUE))>0){
			setgeno<-setgeno[-which(naIdx==TRUE)]
			hmatSample<-hmatSample[-which(naIdx==TRUE)]
		}

		for(j in 1:length(setgeno)){
			cmat<-setgeno[[j]]
			idx.leng0<-which(apply(cmat,2,sum)==0);if(length(idx.leng0)>0){cmat<-cmat[,-idx.leng0]}
			if(ncol(cmat)==1){
				options(warn=-1)
				capture.output(hmat <- homals(data.frame(cmat), eps=1e-05, ndim=1)$obj, file="NUL")
			} else {
				options(warn=-1)
				capture.output(hmat <- homals(data.frame(cmat), eps=1e-05)$obj, file="NUL")
			}
			hmatList[[j]]<-hmat
		}

		names(hmatList)<-names(setgeno)
		comb<-combn(length(setgeno),2)
		ncomb<-ncol(comb)
		p.comb<-rep(NA,ncomb)
	
		for(j in 1:ncomb){
			hmat1<-hmatList[[comb[,j][1]]]; 
			hmat1Sample<-hmatSample[[comb[,j][1]]]
			hmat2<-hmatList[[comb[,j][2]]]; 
			hmat2Sample<-hmatSample[[comb[,j][2]]]
	
			shareSample<-intersect(hmat1Sample, hmat2Sample)
			hmat1<-as.matrix(hmat1[match(shareSample, hmat1Sample),])
			hmat2<-as.matrix(hmat2[match(shareSample, hmat2Sample),])

			can.cor<-try(cancor(hmat1, hmat2)$cor, TRUE)
			if(length(can.cor)<min(ncol(hmat1),ncol(hmat2))){can.cor<-c(can.cor,rep(0,min(ncol(hmat1),ncol(hmat2))-length(can.cor)))}
			capture.output(asym<-p.asym(rho=can.cor, nrow(hmat1), ncol(hmat1),ncol(hmat2),tstat="Wilks"), file="NUL")
			p.comb[j]<-pf(asym$approx,asym$df1,asym$df2,lower.tail=F)[1]	
            

		}
	

		comb1<-apply(as.matrix(comb[,p.comb<0.05/ncomb]),2,function(v)list(na.omit(v)))
		comb1<-lapply(comb1,unlist)
		if(length(which(lapply(comb1,length)==0))>0){comb1<-comb1[-which(lapply(comb1,length)==0)]}

		comb.tab<-table(unlist(comb1))
		comb.replicate<-as.numeric(names(comb.tab)[comb.tab>1])
		if(length(comb.replicate)>0){
			for(j in 1:length(comb.replicate)){
				comb1.idx=which(lapply(comb1,function(v)length(which(v%in%comb.replicate[j]==TRUE)))>0)
				temp.comb1<-comb1[comb1.idx]
				comb1<-comb1[-comb1.idx]
				comb1[[length(comb1)+1]]<-unique(unlist(temp.comb1))
			}
		}
		options(warn=-1)
		unique.comb<-na.omit(unique(unlist(comb1)))
		if(length(unique.comb)>0){	
			setgeno1<-setgeno[-unique.comb]
			hmatSample1<-hmatSample[-unique.comb]
			for(j in 1:length(comb1)){
				intersec<-1:length(y)
				for(k in 1:length(comb1[[j]])){
					intersec<-intersect(intersec,hmatSample[[comb1[[j]][k]]])
				}
				tempgeno<-list()
				for(k in 1:length(comb1[[j]])){
					tempgeno[[k]]<-as.matrix(setgeno[[comb1[[j]][k]]])[match(intersec,hmatSample[[comb1[[j]][k]]]),]
				}					
				setgeno1[[length(setgeno)-length(unique.comb)+j]]<-do.call(cbind,tempgeno)
				hmatSample1[[length(setgeno)-length(unique.comb)+j]]<-intersec
				names(setgeno1)[length(setgeno)-length(unique.comb)+j]<-paste(names(setgeno)[comb1[[j]]],collapse="_")
			}	
		} else {
			setgeno1<-setgeno
			hmatSample1<-hmatSample
		}
	
		setsnp0<-lapply(setgeno1,colnames)
		genesnp0<-unlist(lapply(setsnp0,function(v)paste(v,collapse=",")))
		setmaf0<-lapply(setgeno1,function(v)apply(v,2,get.maf))
		setmac0<-lapply(setgeno1,function(v)apply(v,2,function(v1)length(which(v1>0))))
		setcoef0<-mapply(function(v1,v2)list(get.coef(lm(y~.,data=data.frame(v2)))), hmatSample1, setgeno1)
		setbeta0<-lapply(setcoef0,function(v)v[,1])
		setse0<-lapply(setcoef0,function(v)v[,2])
		setpval0<-lapply(setcoef0,function(v)v[,4])
	

		if(genePerm==TRUE){
			rsltP0<-mclapply(1:length(setgeno1),function(k)QTestPerm.one(y=y,newgeno=as.matrix(setgeno1[[k]]),covadat=NULL, nPerm=nPerm, STT=STT,weight=weight,mc.cores=mc.cores), mc.cores=1)
			rslt.set.q[[i]]<-data.frame(gene=names(setgeno1),size = unlist(lapply(setgeno0,function(v)ncol(as.matrix(v)))),do.call(rbind,rsltP0))
			colnames(rslt.set.q[[i]])<-c("gene","size","p.Q1","p.Q2","p.Q3","p.Q1P","p.Q2P","p.Q3P")
		}else{
			rslt.set.q[[i]]<-data.frame(gene=names(setgeno1),size = unlist(lapply(setgeno1,function(v)ncol(as.matrix(v)))),
			snp=genesnp0, 
			t(mapply(function(v1,v2)c(QTest.one(y=y[v2],newgeno=as.matrix(v1),STT=STT,weight=weight)$p.value, SKAT.one(y[v2],geno=as.matrix(v1)),SKAT.one(y[v2],geno=as.matrix(v1),weights.beta=c(1,1))),setgeno1,hmatSample1)))
			colnames(rslt.set.q[[i]])<-c("gene","size","snp","p.Q1","p.Q2","p.Q3","p.SKAT","p.SKATO","p.SKATc","p.SKATOc")
		}
	
		rslt.gene[[i]]<-data.frame(method="GSQ", gene=rep(names(setgeno1), unlist(lapply(setgeno1,function(v)ncol(as.matrix(v))))), snp=unlist(setsnp0), maf=unlist(setmaf0), mac=unlist(setmac0), beta=unlist(setbeta0) , se=unlist(setse0), pval=unlist(setpval0))
	
		## method=="BQ"
		setsnp1<-lapply(setgeno,colnames)
		setgeno1<-lapply(setgeno,function(v)apply(v,1,sum))
		setgeno2<-do.call(cbind,setgeno1)
		colnames(setgeno2)<-set.list1[[conv.idx]]
		if(length(covadat)!=0){resid<-try(resid(glm(y~.,data=data.frame(covadat),na.action=na.exclude)),TRUE)}
		if(length(covadat)==0){resid<-try(resid(glm(y~1,na.action=na.exclude)),TRUE)}
		S<-setgeno2
		fit<-try(glm(resid~.,data=data.frame(S)),TRUE)
		na.S<-try(which(is.na(coef(fit)[-1])==TRUE),TRUE)
		if(length(na.S)>0){
			S<-try(as.matrix(S[,-na.S]),TRUE)
		setsnp0<-setsnp0[-na.S]
		fit<-try(glm(resid~.,data=data.frame(S)),TRUE)
		}
		coef<-try(coef(summary(fit))[-1,1:2],TRUE)
		if(length(coef)!=2){beta1<-coef[,1];se1<-coef[,2]}
		if(length(coef)==2){beta1<-coef[1];se1<-coef[2]}	
		p2<-pchisq((beta1/se1)^2,df=1,lower.tail=FALSE)
		genesnp1<-unlist(lapply(setsnp1,function(v)paste(v,collapse=",")))
		rslt.set.b[[i]]<-data.frame(gene=colnames(S),size = unlist(lapply(setsnp1,length)), snp=genesnp1, beta=beta1, se=se1, pvalue= p2)

		setmaf1<-lapply(setgeno,function(v)apply(v,2,get.maf))
		setmac1<-lapply(setgeno,function(v)apply(v,2,function(v1)length(which(v1>0))))
		setlm1<-lapply(setgeno,function(v)lm(y~.,data=data.frame(v)))
		setcoef1<-lapply(setlm1,get.coef)	
		setbeta1<-lapply(setcoef1,function(v)v[,1])
		setse1<-lapply(setcoef1,function(v)v[,2])
		setpval1<-lapply(setcoef1,function(v)v[,4])	
		rslt.gene.b<-data.frame(  method="GSB", gene=rep(names(setgeno), unlist(lapply(setgeno,function(v)ncol(as.matrix(v))))), snp=unlist(setsnp1), maf=unlist(setmaf1), mac=unlist(setmac1), beta=unlist(setbeta1) , se=unlist(setse1), pval=unlist(setpval1))
		rslt.gene[[i]]<-rbind(rslt.gene[[i]],rslt.gene.b)
	}

	rslt.set.q1<-data.frame(set=rep(names(idx.accept),unlist(lapply(rslt.set.q,nrow))),do.call(rbind,rslt.set.q))
	rslt.set.b1<-data.frame(set=rep(names(idx.accept),unlist(lapply(rslt.set.b,nrow))),do.call(rbind,rslt.set.b))
	rslt.gene1<-data.frame(set=rep(names(idx.accept),unlist(lapply(rslt.gene,nrow))), do.call(rbind,rslt.gene))
	write.table(as.matrix(rslt.set.q1),paste("withinSet_GSQ_",outname,".rslt",sep=""),quote=F,row.names=F,sep="\t")
	write.table(as.matrix(rslt.set.b1),paste("withinSet_GSB_",outname,".rslt",sep=""),quote=F,row.names=F,sep="\t")
	write.table(as.matrix(rslt.gene1),paste("withinSet_snp_level_",outname,".rslt",sep=""),quote=F,row.names=F,sep="\t")
	return(list(GSQ=rslt.set.q1, GSB=rslt.set.b1, variantLevel=rslt.gene1))
}



get.coef<-function(fit){
	beta0<-coef(fit)[-1]
	coefs<-matrix(coef(summary(fit))[-1,],ncol=4)
	idx.alias<-which(is.na(beta0)==TRUE)	
	if(length(idx.alias)>0){
	coefs1<-se1<-p1<-list()
	i=j=1
	while(i<=length(beta0)){
	if(i%in%idx.alias) {
	coefs1[[i]]<-rep(NA,4)
	i<-i+1
	} else {
	coefs1[[i]]<-coefs[j,]
	i<-i+1
	j<-j+1
	}
	}
	coefs2<-do.call(rbind,coefs1)
	rownames(coefs2)<-names(beta0)
	} else {
	coefs2<-coefs
	}
	return(coefs2)
}




