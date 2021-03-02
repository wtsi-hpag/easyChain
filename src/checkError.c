#include "sys/sysinfo.h"
#include"sys/time.h"
#include<dirent.h>
#include <sys/mman.h> 
#include<sys/types.h>
#include<sys/stat.h>
#include<stdlib.h> 
#include<memory.h> 
#include<iostream> 
#include<string.h>
#include <fcntl.h>
#include<fstream>
#include <map>
#include <vector>
#include <pthread.h>
#include <unistd.h>
#include <error.h>
#define _MAX_INT_DIG 32
#define CHR_NUM 24
#define BIN_FILE_NUM 24
#define _GLIBCXX_USE_CXX11_ABI 0   // close c++11

using namespace  std;

string sHelp = "Usage: ./checkError [VCF_Files_folder] [reference_VCF_Files_folder] [chain_file] [output_folder]\n";
string sExample = "Example: ./checkError grch37_vcf grch38_vcf hg19ToHg38.over.chain result \n";
string sHelpInfo = "[Usage]: ./checkError [VCF_Files_folder] [reference_VCF_Files_folder] [chain_file] [output_folder] \n\
[Example]: ./checkError grch37_vcf grch38_vcf hg19ToHg38.over.chain output_result \n\n\
   grch37_vcf            - The folder which contains some or all the VCF files for GRCh37\n\
   grch38_vcf            - The folder which contains ALL the VCF files for GRCh38\n\
   hg19ToHg38.over.chain - The chain file selected\n\
   output_result         - The folder with the output results \n\n\
[Note]:  Five files will be generated after processing for each input VCF file, they are:\n\
 [1]  xxx_SNP.bed:                   - The bed file extracted from xxx.vcf with tag \"VT=SNP\".\n\
 [2]  xxx_SNP_genegos.bed:           - The file after coordinate conversion.\n\
 [3]  xxx_SNP_genegos.unmap:         - The content that could not be converted.\n\
 [4]  xxx_SNP_genegos_error.dat:     - The file contains all the error sites.\n\
 [5]  xxx_SNP_genegos_error_db.txt:  - The file contains all error sites in import databases.\n";

typedef  unsigned	int  Uint;

struct ChainBag
//class ChainBag
{
	public:
	bool bReverse;
	//bool bHasSun;
	int  nSrcStart;
	int  nSrcEnd;
	int  nDesStart;
	int  nDesEnd;
	int  nSpanChr; 
	string strDesChr;
	//nSpanChr
	ChainBag()
	{
		bReverse = false;
		//bHasSun = false;
		nSrcStart = -1;
		nSrcEnd = -1;
		nDesStart =-1;
		nDesEnd = -1;
		nSpanChr = 0;
		strDesChr = "";
	}
	~ChainBag()
	{
		bReverse = false;
		//bHasSun = false;
		nSrcStart = -1;
		nSrcEnd = -1;
		nDesStart =-1;
		nDesEnd = -1;
		nSpanChr = 0;
		strDesChr = "";
	}
};

struct  CmySection
//class CmySection
{
	public:
	std::map<int,int> map_Chain_1V1;
	std::map<int,ChainBag> map_Chain;
	CmySection()
	{
		map_Chain_1V1.clear();
		map_Chain.clear();
	}
	~CmySection()
	{
		map_Chain_1V1.clear();
		map_Chain.clear();
	}
};

class THREAD_FILEMMAP_ARG
{
public:
	FILE *p_snp;
	FILE *p_indel;
	char *pStart;
	long nlen;
	string sResult; 
	string sFail;
	CmySection *makeUp;
	map<string,CmySection> *map_Section; // convert via stand chain
	
	THREAD_FILEMMAP_ARG()
	{
		nlen = -1;
	  sResult = "";
    sFail  = "";  
    pStart = NULL;
	  makeUp = NULL;
	  map_Section = NULL;
	}
	~THREAD_FILEMMAP_ARG()
	{
		sResult = "";
    sFail  = ""; 
    pStart = NULL;
    makeUp = NULL;
	}
};


struct THREADARG 
{
	int nNum ;
	map<Uint,Uint> *Map_PosBag;
	string sFileName;
	THREADARG()
	{
		nNum = 0;
		sFileName = "";
		Map_PosBag =NULL;
	}
};

struct THREADARGESDB
{
	multimap<string,string> *pDAT_V2;
	string sFile;
	bool bInit;
	THREADARGESDB()
	{
		sFile ="";
		bInit = false;
	}
};

 THREADARGESDB Thrdbag[CHR_NUM];
 pthread_t Init_thread[CHR_NUM];              //MT initialize
 multimap<string,string> Map_DAT_V2[CHR_NUM]; //data base GRCh38 of clinvar gwas hgmd omim



typedef std::map <int,int>::iterator Map_IntInt_iterator;
typedef std::map <int,ChainBag>::iterator Map_Chain_iterator;
typedef  pair<map<int,ChainBag>::iterator,bool> MapRet;
typedef std::map<string,CmySection>::iterator mapChr_Chains_iterator;
typedef std::map<string,std::map<int,ChainBag> >::iterator ChrChain_iterator;




/*****************************************************************************
** Function:     int2string(int _Val)
** Create Date:  2018.7.17
** Modify Time:  2018.7.17
** Author:        LYJ
** Version:      1.0
*******************************************************************************/

string int2string(int _Val)
{	
	char _Buf[2 * _MAX_INT_DIG];
  snprintf(_Buf, sizeof(_Buf), "%d", _Val);
	return (string(_Buf));
}

int str2int(string str)
{
	int num = atoi( str.c_str() );
	return num;
}

Uint stoUint(string str)
{
	unsigned int result(0);  //4294967296（=2^32-1）
	if(str.length() >10)
		return 0;
	for (int i = str.size()-1;i >= 0;i--)
	{
		unsigned int temp(0),k = str.size() - i - 1;
		if (isdigit(str[i]))
		{
			temp = str[i] - '0';
			while (k--)
				temp *= 10;
			result += temp;
		}
		else
			break;
	}
	return result;
}


string Uint2string(Uint _Val)
{	
	char _Buf[2 * _MAX_INT_DIG];
  snprintf(_Buf, sizeof(_Buf), "%u", _Val);
	return (string(_Buf));
}



/*****************************************************************************
** Function:     is_file_exist(string sFilename)
** Create Date:  2018.7.17
** Modify Time:  2018.7.17
** Author:        LYJ
** Version:       1.0
*******************************************************************************/

 bool is_file_exist(string sFilename)
 {
 	if( 0 != access(sFilename.c_str(),F_OK))
 		return false;
  return true;
 }
 
 
bool is_dir_exist(string sDir)
 {
 	
 	  DIR *dirptr=opendir(sDir.c_str());
 	  
 	  if( NULL == dirptr)
 	  	return false;
 	  
 	  closedir(dirptr);
 	  return true;
 }
 
 
int CreateDir(const char *sPathName)  
{  
  char   DirName[512];  
  strcpy(DirName,sPathName);  
  int   i,len   =   strlen(DirName);  
  if(DirName[len-1]!='/')  
  strcat(DirName,"/");  
  len  =   strlen(DirName);  
   
  for(i=1; i<len; i++)  
  {  
    if(DirName[i]=='/')  
    {  
       DirName[i] = 0;  
      if( access(DirName,F_OK)!=0   )  
       {  
        if(mkdir(DirName,0755)==-1)  
         {   
          perror("mkdir error");   
          return   -1;   
         }  
       }  
     DirName[i]   =   '/';  
    }  
  }  
   return   0;  
} 
 
 
 /*****************************************************************************
** Function:      bWritResult
** write str to file
** Create Date:  2018.7.04
** Modify Time:  2018.7.04
** Author:       LYJ
** Version:      1.0
*******************************************************************************/
bool bWritResult(FILE * &pFile,string *pstrInfo)
{
	size_t nLen = pstrInfo->length();
	size_t nWriten= 0;
	nWriten = fwrite (pstrInfo->c_str(), 1, nLen, pFile);

	if (nLen != nWriten)
	{
		printf("write file length error src_len: %ld  written_len: %ld \n",nLen ,nWriten);
		return false;
	}
	return true;
}

/*****************************************************************************
** Function:      Replace_char(string& str,const string&old_value,const string& new_value)
** Create Date:  2018.7.05
** Modify Time:  2018.7.05
** Author:        LYJ
** Version:       1.0
*******************************************************************************/
int Replace_char(string &str,const char &old_value,const string& new_value)
{
	int nRe = 0;
	for(int pos=0;pos!=-1;pos+=new_value.length())
	{
		if((pos=str.find(old_value,pos))!=-1)
		{
			//str.replace(pos,old_value.length(),new_value);
			str.replace(pos,1,new_value);
			nRe++;
		}
		else 
			break;
	}
	return nRe;
}

int Replace_char(string &str,const string &old_value,const string& new_value)
{
	int nRe = 0;
	for(int pos=0;pos!=-1;pos+=new_value.length())
	{
		if((pos=str.find(old_value,pos))!=-1)
		{
			str.replace(pos,old_value.length(),new_value);
			nRe++;
		}
		else 
			break;
	}
	return nRe;
}

bool  getCurrentPath(string &sPath )
{
	char current_path[1024]; 
  int cnt = readlink("/proc/self/exe",current_path, 1024); 
  if (cnt < 0 || cnt >= 1024) 
  { 
    printf("***Error***\n"); 
     return false;
  } 
  int i; 
  for (i = cnt; i >=0; --i) 
  { 
    if (current_path[i] == '/') 
    { 
        current_path[i] = '\0'; 
        break; 
    }
  } 
  sPath = current_path ;
  return true;
}


long GetFileSize(string sfname)
{
	 int fd;  
	 struct stat sb;  
	
	 if ((fd = open(sfname.c_str(), O_RDONLY)) < 0)
	 	return -1;
	  
	 if ((fstat(fd, &sb)) == -1)
	  	return -2;
	  	
	 long  nFileSize  = sb.st_size;
	 close(fd);
	  
	return nFileSize;
}

long GetFileSize(int fd)
{ 
	 struct stat sb;  
	 if ((fstat(fd, &sb)) == -1)
	  	return -1;
	return sb.st_size;
}

 static int GetChr(string strch)
  {
  	Replace_char(strch,"chr","");
  	if(strch == "" || strch.find_first_not_of("0123456789xyXY")!= -1)
  		return -1;
  		
		int nNum = 0;
		if (strch == "X"||strch == "x")
			return 23;
		else if (strch == "Y"||strch == "y")
			return 24;
		else
			nNum = str2int(strch);
	  
	  return nNum;
  }
  static string GetChr(int nChr)
  {
  	string strch = "chr";
  	
  	if(nChr == 23)
  		strch = "chrX";
  	else if(nChr == 24)
  			strch = "chrY";
  	else
  		  strch += int2string(nChr);
	  return strch;
  }
  

  static int GetPos(string strPos)
  {
  	if(strPos == "" || strPos.find_first_not_of("0123456789")!= -1)
  		return -1;
		return str2int(strPos);
  }
  
int  getAllFiles(string path, vector<string>& files) 
{
	int nRe = 0;
	DIR *dirp;
  struct dirent *dp;
  dirp = opendir(path.c_str());
  while ((dp = readdir(dirp)) != NULL) 
    {
    	  if(strcmp(dp->d_name,".")==0 || strcmp(dp->d_name,"..")==0)    ///current dir OR parrent dir
    	  	  continue;
        
        string strName(dp->d_name);
        strName = path + "/" + strName;
        Replace_char(strName,"//","/");
        files.push_back(strName);
        nRe++ ;
    }
    closedir(dirp);
	return nRe;
}
 
 
  
 
/*****************************************************************************
** Function:    GetPosSectionMakeUp(int nSpos,CmySection makeUp)
** Get coordinat conversion pos
** Create Date:  2018.7.17
** Modify Time:  2018.7.17
** Author:        LYJ
** Version:      1.0
*******************************************************************************/
static inline int GetPosSectionMakeUp(int nSpos,CmySection *pmakeUp)
{
	int nRe = -1;
	Map_IntInt_iterator itInt = pmakeUp->map_Chain_1V1.find(nSpos);
	if (itInt != pmakeUp->map_Chain_1V1.end())
	{
		nRe = itInt->second;
		return nRe;
	}
	
	Map_Chain_iterator it = pmakeUp->map_Chain.begin();
	for (;it != pmakeUp->map_Chain.end();it++)
	{
		if (nSpos >=it->first && nSpos<= it->second.nSrcEnd)
		{
			bool bre= it->second.bReverse;
			int nDif = nSpos - it->first;
			if (!bre)//"+"
				nRe = it->second.nDesStart + nDif;
			else
				nRe = it->second.nDesStart - nDif; 
			break;
		}
		else if (nSpos <it->first && nSpos< it->second.nSrcEnd)  
		{
			break;
		}

	}
	return nRe;

}

static inline int GetPosSectionMakeUp(int nSpos,CmySection *pmakeUp,string &sDesChr)
{
	int nRe = -1;
	if(nSpos<0)
		return nRe;
		
	Map_IntInt_iterator itInt = pmakeUp->map_Chain_1V1.find(nSpos);
	if (itInt != pmakeUp->map_Chain_1V1.end())
	{
		nRe = itInt->second;
		return nRe;
	}
	//std::map<int,ChainBag> &map_Chain = pmakeUp->map_Chain;
	Map_Chain_iterator it = pmakeUp->map_Chain.begin();
	for (;it != pmakeUp->map_Chain.end();it++)
	{
		if (nSpos >=it->first && nSpos<= it->second.nSrcEnd)
		{
			if(""!= it->second.strDesChr)
			  sDesChr = it->second.strDesChr;
			
			if (!it->second.bReverse)//"+"
				nRe = it->second.nDesStart + nSpos - it->first;
			else
				nRe = it->second.nDesStart - nSpos + it->first; 
			break;
		}
		else if (nSpos <it->first && nSpos< it->second.nSrcEnd)
		{
			break;
		}

	}
	return nRe;

}

///////////////////////modify by lyj 11.22 to insure stand chains 
static inline int GetPosSectionMakeUp(int nSpos,CmySection *pmakeUp,string &sDesChr,bool &bRevs)
{
	int nRe = -1;
	if(nSpos<0)
		return nRe;
		
	Map_IntInt_iterator itInt = pmakeUp->map_Chain_1V1.find(nSpos);
	if (itInt != pmakeUp->map_Chain_1V1.end())
	{
		nRe = itInt->second;
		return nRe;
	}
	//std::map<int,ChainBag> &map_Chain = pmakeUp->map_Chain;
	Map_Chain_iterator it = pmakeUp->map_Chain.begin();
	for (;it != pmakeUp->map_Chain.end();it++)
	{
		if (nSpos >it->first && nSpos<= it->second.nSrcEnd)
		{
			if(""!= it->second.strDesChr)
			  sDesChr = it->second.strDesChr;
			
			bRevs = it->second.bReverse;
			if (!bRevs)//"+"
				nRe = it->second.nDesStart + nSpos - it->first;
			else
				nRe = it->second.nDesStart - nSpos + it->first; 
			break;
		}
		else if (nSpos <it->first && nSpos< it->second.nSrcEnd)
		{
			break;
		}

	}
	return nRe;

}

static inline int GetPosSection(int nSpos,CmySection *pmakeUp)
{
	int nRe = -1;
	if(nSpos<0)
		return nRe;
	
	Map_Chain_iterator it = pmakeUp->map_Chain.begin();
	for (;it != pmakeUp->map_Chain.end();it++)
	{
		if (nSpos >it->first && nSpos<= it->second.nSrcEnd)
		{ 
			nRe  = nSpos - it->first;
			break;
		}
		else if (nSpos <it->first && nSpos< it->second.nSrcEnd)
			break;

	}
	return nRe;
}


/*****************************************************************************
** Function:   bCoordinate(string &sChr,string &sStart,CmySection *pmakeUp)
** convert file via stand chains file 
** Create Date:  2019.12.15
** Modify Time:  
** Author:        LYJ
** Version:       1.0
*******************************************************************************/

  //ordinary VCF mode
 bool inline bCoordinate(string &sChr,string &sStart,CmySection *pmakeUp)
 {
 	   int nStart = GetPos(sStart);
	   if(nStart == -1) {return false;}
	   	
     string sDesChr = "";
     bool bRe_E = false;
     int nDesPos= GetPosSectionMakeUp(nStart,pmakeUp,sDesChr,bRe_E);
     if(nDesPos == -1)
     	return false;
     
     if("" != sDesChr)
     	sChr = sDesChr;
     
     sStart = int2string(nDesPos);
    return true;
 	
 }
 
 //ordinary BED mode
 bool inline bCoordinate(string &sChr,string &sStart,string &sDend,CmySection *pmakeUp)
 {
 	   int nt_start = GetPos(sStart);
	   int nt_end = GetPos(sDend);
	   if(nt_start == -1 || nt_end == -1 ) {return false;}
	   	
	   int nGap = nt_end - nt_start;
     bool bRe_S = false;
     bool bRe_E = false;
     string DesChr1 = "";
     string DesChr2 = "";
     
     int nPosEnd = GetPosSectionMakeUp(nt_end,pmakeUp,DesChr2,bRe_E);
     if(nPosEnd == -1)
     	return false;
     	
     int nPosStart = GetPosSectionMakeUp(nt_start,pmakeUp,DesChr1,bRe_S);
    
     if(nPosStart == -1 || DesChr1 != DesChr2)
     	{
     		if(DesChr2 != "") {sChr = DesChr2;}
     		if(bRe_E)
     			{
     				nPosStart = nPosEnd;
     				nPosEnd += nGap;
     			}
     		else
     			nPosStart = nPosEnd - nGap;
     	}
     	else //nPosStart != -1
     	{
     		 if(DesChr2 != "") {sChr = DesChr2;}
     		 //////////////////////////////////////////////////
     		 if(nGap != abs(nPosStart-nPosEnd)) //错开区间
     		 	{
     		 		 int nGapReal = GetPosSection(nt_start,pmakeUp) +1;
     		 		 if( nGapReal>=nGap )
     		 		 	{
     		 		 		if(bRe_E)
     		 		 			nPosStart = nPosEnd + nGap;
     		 		 		else
     		 		 			nPosStart = nPosEnd - nGap;
     		 		 	}
     		 	}
     		 	//////////////////////////////////////////////////
     		 
     		 if(nPosStart > nPosEnd ) // 落到同一条染色体上
     		 {
     		 			nt_start = nPosStart;
     				  nPosStart = nPosEnd;
     				  nPosEnd = nt_start;
     		 }
     	}
     	sStart = int2string(nPosStart);
     	sDend  = int2string(nPosEnd);
    
    return true;
 	
 }
 
 static inline int nSplitStr2List(string &str,vector<string> &sL,const string splt)
{
	int npS = 0;
	int npE = 0;
	int nRe = 0;
	string strGet = "";
	sL.clear();
	for (;npE != -1;)
	{
		npE = str.find(splt,npS);
		if (npE == -1) // 最后的字串
		{
			strGet = str.substr(npS);//最后
			sL.push_back(strGet);
			nRe++ ;
			continue;
		}
		strGet = str.substr(npS,npE-npS);
		npS = npE +1;
		sL.push_back(strGet);
		nRe ++;
	}
	return nRe;
}
 


 bool bIintChains(string sChainFile,map<string,CmySection> *mapChr_Chains)
{
	ifstream ifs;
	ifs.open(sChainFile.c_str(),ios_base::in);
	if (!ifs)
	{
		cout<<"Read File Error, please check file name is right: "<<sChainFile<<endl;
		return false;
	}

  int ntStart,ntEnd,nqStart,nqEed,nSrcStart,nDesStart,nSrcEnd,nDesEnd,nqSize,nSpan;
  ntStart=ntEnd=nqStart=nqEed=nSrcStart=nDesStart=nSrcEnd=nDesEnd=nqSize=nSpan=0;
	
	string str_Line,sDesChrName,sLastChr,sStand;
	str_Line=sDesChrName=sLastChr=sStand="";
	ChainBag c_bag;
	CmySection sC;
	mapChr_Chains_iterator ChainsIter;
	vector<string> subDis;
	int nline = 0;

	while(!ifs.eof())  
	{  
		getline(ifs,str_Line);
		nline ++;
		if (str_Line == "\r" ||str_Line == "")
			continue;
		
		if (str_Line.find("chain") != -1) // head
		 {
		 	//cout << str_Line << endl;
			if( 0 == nSplitStr2List(str_Line,subDis," ")) //splited by 'space'
				{
					ifs.close();
					cout << "Chain file error line " << nline << endl;
	        return false;
				}
			
			if(sLastChr == "")
				{
					sLastChr = subDis[2];
					mapChr_Chains->insert(std::make_pair(subDis[2],sC));
					//cout << subDis[2] << endl;
					ChainsIter = mapChr_Chains->find(subDis[2]);
				} 
			else if(sLastChr != subDis[2]) // change chain
				{
					ChainsIter = mapChr_Chains->find(subDis[2]);
					if(ChainsIter == mapChr_Chains->end())
					 {
					  	mapChr_Chains->insert(std::make_pair(subDis[2],sC));
					  	ChainsIter = mapChr_Chains->find(subDis[2]);
					  	//cout << subDis[2] << endl;
					  	if(ChainsIter == mapChr_Chains->end())
					  	{
					  		 ifs.close();
					  		 cout << "Chain file error line " << nline << endl;
	               return false;
					  	}
					 }
				}

			nSrcStart = str2int(subDis[5]);
			nDesStart = str2int(subDis[10]);
			//nDesEnd = str2int(subDis[11]);
			nqSize = str2int(subDis[8]);				
			sStand = subDis[9];
			c_bag.strDesChr = (subDis[2] == subDis[7])?(""):(subDis[7]);
			c_bag.bReverse = (sStand == "-")?(true):(false);
		 }
		else // normal three or one    // "167417 50000 80249" or "40302"
			{
				Replace_char(str_Line," ","\t");// ensemble is space ; ucsc is \t
				int nCell = nSplitStr2List(str_Line,subDis,"\t");
				nSpan = str2int(subDis[0]);
			  c_bag.nSrcStart = nSrcStart;
			  c_bag.nSrcEnd = nSrcStart + nSpan;
			  nDesEnd = nDesStart + nSpan;
			  if (c_bag.bReverse) //"-"
			   {
				    c_bag.nDesStart = nqSize - (nDesStart);
				    c_bag.nDesEnd = nqSize - nDesEnd ;
			   }
			 else
			   	{
				    c_bag.nDesStart = nDesStart;
				    c_bag.nDesEnd = nDesEnd;
			    }
			  if( 3 == nCell )
			  	{ 
			        nSrcStart = c_bag.nSrcEnd + str2int(subDis[1]);
			        nDesStart = nDesEnd + str2int(subDis[2]);
			    }
			    MapRet ret =  ChainsIter->second.map_Chain.insert(std::make_pair(c_bag.nSrcStart,c_bag));
			    	
			    //if(c_bag.bReverse == false)
			   // 	cout << c_bag.nSrcStart << "\t" << c_bag.nSrcEnd << "\t+\t" << c_bag.nDesStart << "\t" << c_bag.nDesEnd << endl;
			   // else
			   // 	cout << c_bag.nSrcStart << "\t" << c_bag.nSrcEnd << "\t-\t" << c_bag.nDesStart << "\t" << c_bag.nDesEnd << endl;
			    	
				  if (!ret.second)
				   {
					    //ifs.close();
					    //cout <<sLastChr << endl;
					    //cout << c_bag.nSrcStart << "\t" << c_bag.nSrcEnd << "----->\t" << c_bag.nDesStart << "\t" << c_bag.nDesEnd << endl;
	            //return false;
				   }
			}
	}
	
	ifs.close();
	return true;
}


  
  void * TCnvtBed_stand(void* pT)
	{
		  THREAD_FILEMMAP_ARG *pBag = (THREAD_FILEMMAP_ARG *)pT;
			register char *pChar = pBag->pStart;
	    register long nlen = pBag->nlen;
	    map<string,CmySection>* pChr_Sec = pBag->map_Section;
	    mapChr_Chains_iterator it;
	    int nChr = -1;
	    int nspanChr = -1;
	    int nPos = -1;
	    int nTCout = 0;
	    register int  nslen =0;
	    register int  nLine_len =0;
	    int  nreflen =0;  // 标记chr + pos 长度；
	    string sChr = "";
	    string sPos = "";
	    string sPosStart = ""; // posend
	    string sPosEnd = ""; // posend
	    string sLine = "";
	    string sLastChar = "";
	    bool  bConvert = false;
	    int nPosStart = -1;
	    int nPosEnd = -1;
	    string DesChr1 = "";
	    string DesChr2 = "";
	    
	    
	    pBag->sResult.reserve(nlen*1.1); // 设置好文件内存使用大小
	    char *pS = pChar;
	    char *pStart_Line = pChar;
	    char *pStart_REF = pChar;
	    
	   for(;nlen >=0 ;nlen--)
	   //for(;nlen >=0 ;nlen--)
	    {
	    	switch(*pChar++)
	    	{
	    		case'\n':
	    			{
	    				if(bConvert)
	    					{
	    						sLine.assign(pStart_REF, nLine_len-nreflen+1);
	    						//cout << sLine << endl;
	    						sLine = sChr + "\t" + sPosStart + "\t"  + sPosEnd + "\t" + sLine;
	    						pBag->sResult += sLine;
	    					}
	    				 else
	    				 	{
	    				 		sLine.assign(pStart_Line, nLine_len+1);
	    		        //if(sLine.find("#")== 0)
	    		        if(sLine[0]== '#')
	    		        	pBag->sResult += sLine;
	    		        else
	    		        	pBag->sFail += sLine;
	    				 	}
	            nTCout = nslen = nreflen = nLine_len = 0;
	            sChr = sPos = "";
	            pS = pStart_Line = pChar;
	            bConvert = false;
	    				break;
	    			}
	    		case '\t':
	    			{
	    				nTCout ++;
	    				nLine_len++;
	    				
	    				if(nTCout<=3)
	    				{
	    			      switch(nTCout)
	    			    	{
	    					     case 1:
	    						    sChr.assign(pS,nslen);
	    					      break;
	    					     case 2:
	    						     {
	    							     sPosStart.assign(pS,nslen);
	    							     if( sChr != sLastChar )
	    								     {
	    									     sLastChar = sChr;
	    									     it = pChr_Sec->find(sChr);
	    							         if(it == pChr_Sec->end())// chrx //x
	    								          it = pChr_Sec->find("chr"+sChr);
	    								       if(it == pChr_Sec->end())// x
	    								  	    {
	    								  		     string s = sChr;
	    								  		     Replace_char(s,"chr","");
	    								  		     it = pChr_Sec->find(s);
	    								  	    }    
	    								     }	
	    						       break;
	    						     } 
	    					      case 3:
	    						    {
	    							     sPosEnd.assign(pS,nslen);
	    							     //cout << sChr << "\t" << sPosStart << "\t" << sPosEnd << endl;
	    							   if(it != pChr_Sec->end()){
	    								    if(bCoordinate(sChr,sPosStart,sPosEnd,&it->second))
	    							      {
			                      pStart_REF = pChar;
			                      nreflen = nLine_len;
			                      bConvert = true;
			                  }
	    								
	    							  }
	    						   break;
	    						    }
	    					    default:
	    							   break;
	    						
	    				    }
	    						
	    				}
	    				pS = pChar;
	    				nslen =0;
	    				break;
	    			}
	    	  default:
	    	  	{
	    	  		nslen ++;
	    	  		nLine_len++;
	    	  		break;
	    	  	}
	    	  	
	    	} ///////////end switch 
	    }////////////end for
	    
	   ////////////////////////最后一行没有\n的数据////////////////////////////// 
	   if(nLine_len >0)
	   	{
	   	  if(bConvert)
	   	  	{
	    	   sLine.assign(pStart_REF,nLine_len-nreflen-1);
	    	   sLine = sChr + "\t" + sPosStart + "\t"  + sPosEnd + "\t" + sLine;
	    	   pBag->sResult += sLine;
	        }
	      else
	      	{
	    	    sLine.assign(pStart_Line, nLine_len-1);
	         if(sLine.find("#")== 0)
	    		   pBag->sResult += sLine;
	    	   else
	    		   pBag->sFail += sLine;
	        }
	   	}
	   //////////////////////////////////////////////////////// 
	    
	   string strRe = " finish \n";
	   char *cSt =(char *) strRe.data();
	  return (void *)cSt;
  }
  


////////////////////convert file via stand chain//////////////////////////////////////////////////////////////////
bool Cnvt_File_Stand(string sSrcfile,string sChainFile,string sDespath,string &sError)
{
	bool bRe = true;
	int fd;  
	char *mapped;
	bool bVcf = true;
	long nFilesize =0;
	FILE *pFile;
	FILE *pFSpecial;
	string strSection = "";
	sError = "";
	
	int ncpuNum = get_nprocs();
	ncpuNum = (ncpuNum>1)?(ncpuNum-1):1;
	
			
	pthread_t Init_thread[ncpuNum];
	THREAD_FILEMMAP_ARG  Thread_Check[ncpuNum];
	     
	 if ((fd = open(sSrcfile.c_str(), O_RDONLY)) < 0)
	  {  
	      perror("open"); 
	      sError = "open file error: " + sSrcfile;
	      return false;  
	  } 
	  
	  nFilesize = GetFileSize(fd);
	  if(nFilesize < 5)
	  	{  
	  		close(fd);
	  		cout << sSrcfile << "is empty, please check\n";
	  		return false;
	  	}
	  	
	  if ((mapped = (char*)mmap(NULL, nFilesize, PROT_READ, MAP_SHARED, fd, 0)) == (void *)-1) 
	   {  
	   	    close(fd);
	        perror("mmap"); 
	        sError = "mmap file error: " + sSrcfile;
	        return false; 
	   }
	  close(fd);
	///////////////////////////////////////////intital convert section via stand chains/////////////////////////////////////
	 map<string,CmySection> mapChr_Chains;
	 if( false == bIintChains(sChainFile,&mapChr_Chains))
	 	{
	 		    close(fd);
	        sError = "initialize chains error\n"; 
	        return false; 
	 	}
	 	
	 	//cout << "begin change \n";
	 	//map<string,CmySection>::iterator iiit = mapChr_Chains.begin();
	 	
	//////////////////////////////////////////////MAKE DES FILES////////////////////////////////////////////////////////////
	
	string sFilename = sSrcfile.substr(sSrcfile.rfind('/')+1,sSrcfile.length() - sSrcfile.rfind('/') -1 );
	string mapFile,postfix;
	mapFile = sFilename.substr(0,sFilename.rfind('.'));
	postfix = sFilename.substr(sFilename.rfind('.'));
	mapFile += "_genegos" + postfix;
	mapFile = sDespath + "/" +  mapFile;
	Replace_char(mapFile,"//","/");

	
	if ((pFile = fopen(mapFile.c_str(),"w+"))==NULL)
	{
		sError = "Creat file error: " + mapFile;
		return false;
	}
	sError = mapFile;
	cout << "After coordinate conversion bed file: " << mapFile << "\n";
	
	string smark_file = sFilename.substr(0,sFilename.rfind('.'));
	smark_file += "_genegos.unmap";
	smark_file = sDespath + "/" +  smark_file;
	Replace_char(smark_file,"//","/");
	
	if ((pFSpecial = fopen(smark_file.c_str(),"w+"))==NULL)
	{
		sError = "Creat file error: " + smark_file ;
		if (fclose (pFile) != 0)
			perror("Error occurs when close file"); 
		
		return false;
	}
	
	 long nPerThread = (nFilesize >40960)? (nFilesize/ncpuNum):nFilesize; //split mem
	 long nLeft = nFilesize;
	 char *pStart = mapped;
	 long  nlen = 0;
	 int nRealThread = 0;
	 
	  for(int i= 0;i< ncpuNum && nLeft >0;i++)
	  {
      nRealThread ++;
			Thread_Check[i].pStart = pStart;
			Thread_Check[i].makeUp  = NULL ;
			Thread_Check[i].map_Section = &mapChr_Chains;
			if(nLeft <= nPerThread)
			{
				Thread_Check[i].nlen = nLeft;
				nLeft = 0;
				break; 
			}
				  
			nLeft -= nPerThread;
		  pStart += nPerThread;
		  nlen = nPerThread;
			for(;(*pStart) != '\n'&&nLeft>0;pStart++)
			{
				nLeft --;
				nlen ++;
			}
			
			 pStart++;
		 	 nLeft --;
		 	 nlen++;
			 Thread_Check[i].nlen = nlen;
	  }
	  
	  
	  //nRealThread = 2;
	   for (int i=0;i< nRealThread;i++)
		 {
			int nerror = 0;
			nerror = pthread_create(&Init_thread[i],NULL,TCnvtBed_stand,(void *)&Thread_Check[i]);
			if (0!=nerror)
			{
				sError = "Creat thread error \n";
	      if (fclose (pFile) != 0)
	      	perror("Error occurs when close file");
			  return false;
			}
		 } 
		 
		for (int i = 0; i <nRealThread ;i++)
		{
				char *ps = NULL;
				int nerror = pthread_join(Init_thread[i],(void **)&ps);
				if (0!=nerror)
				{
					sError = "Wait thread error \n";
					printf("pthread_join  %s \n", strerror(nerror));
					return false;
				}
				else 
				{
					bWritResult(pFile,&Thread_Check[i].sResult);
					bWritResult(pFSpecial,&Thread_Check[i].sFail);
				}
		}
		
	  if ((munmap((void *)mapped, nFilesize)) == -1)
	  	{
	  		perror("munmap"); 
	  		sError = "unmap file error: " ;
	  		bRe = false;
	  	}
	      

	  if (fclose (pFile) != 0 )
	  	{
	  		perror("Error occurs when close file"); 
	  		sError = "Error occurs when close file: " +  mapFile;
	  	}
	  	
	  if (fclose (pFSpecial) != 0)
	  	{
	  		perror("Error occurs when close file"); 
	  		sError = "Error occurs when close file: " +  smark_file;
	  	}
	  	
	return bRe;
}


 static int n_PickInfosnp(string &str)
{
	bool bNormal = false;
	bool bSNP = true;
	int nposEnd = 0;
	int npos1 = 0;// 自串头位置
	int npos2 = 0;// 自串尾巴位置
	int nReturn = 0;// 返回值 0:snp  1:indel -1:error
	
	npos2 = str.find("	",0);
	string strGet;
	if (npos2!= -1)
	{
		/********************************CHORM**********************************/
			npos1 = npos2+1;
		/********************************POS*************************************/

			npos2 = str.find("	",npos1);
			string sPos= str.substr(npos1,npos2-npos1);
			int  nPos = str2int(sPos)+1;
			sPos = int2string(nPos);
			//strGet += sPos;
			strGet = "chr" + str.substr(0,npos2+1) + sPos + "\t";
			npos1= npos2+1;
		/********************************ID*************************************/
			npos2 = str.find("	",npos1);
	  	string sID= str.substr(npos1,npos2-npos1);
			strGet += sID;
			//strGet = str.substr(npos1,npos2-npos1);
			npos1= npos2+1;
		/********************************REF***************************************/
			npos2 = str.find("	",npos1);
			//nReflen = npos2 - npos1;
			//string sREF= str.substr(npos1,npos2-npos1);
			//if (npos2-npos1>=2)
			//	bNormal = true;
			//strGet += sREF;
			//strGet = str.substr(npos1,npos2-npos1);
			npos1= npos2+1;
		/********************************ATL***************************************/
			npos2 = str.find("	",npos1);
			//sATL= str.substr(npos1,npos2-npos1);
			/*if (npos2-npos1>=2)
			bNormal = true;*/
			//strGet += sATL;
			//strGet = str.substr(npos1,npos2-npos1);
			npos1= npos2+1;
			nposEnd = npos1;// 不要pass
		/********************************QUAL***************************************/
			npos2 = str.find("	",npos1);
			//string sQUAL= str.substr(npos1,npos2-npos1);
			//strGet += sQUAL;
			//strGet = str.substr(npos1,npos2-npos1);
			npos1= npos2+1;
		/********************************FILTER***************************************/
			npos2 = str.find("	",npos1);
		//	string sFILTER= str.substr(npos1,npos2-npos1);
			if (npos2-npos1 != 4)
				bNormal = true;
			//strGet += sFILTER;
			//nposEnd = npos2 +1;   // 要pass字段
			npos1= npos2+1;
			//if ( "PASS" != sFILTER  )
			//cout << sFILTER <<endl;

		if (str.find("VT=SNP",npos1) == -1)//
		{ 
			if(sID.find("rs") == -1)
				nReturn = -1;
			else
				nReturn = 1;
		}
		
		str = strGet + "\n";
		}
	else 
	{
		return -1;
	}
	return nReturn;
}


 static int n_PickRS(string &str,Uint &nPos,Uint &urs,int &nChr)
{
	int npos1 = 0;// 自串头位置
	int npos2 = 0;// 自串尾巴位置
	int nReturn = 0;// 返回值 0:snp  1:indel -1:error
	
	npos2 = str.find("\t",0);
	if (npos2!= -1)
	{
		/********************************CHORM**********************************/
		  string sChr= str.substr(0,npos2);
		  Replace_char(sChr,"chr","");
		  nChr = str2int(sChr);
			npos1 = npos2+1;
		/********************************POS*************************************/
			npos2 = str.find("\t",npos1);
		  string sPos= str.substr(npos1,npos2-npos1);
		  nPos = stoUint(sPos);
			npos1= npos2+1;
		/********************************ID*************************************/
			npos2 = str.find("\t",npos1);
	  	string sID= str.substr(npos1,npos2-npos1);
	  	if(sID.find_first_not_of("rs0123456789") != -1 )
	  		return -1;
	  		
	  	if(Replace_char(sID,"rs","")!= 1)
	  		return -1;
	  	
	  	urs = stoUint(sID);
      /* if (str.find("VT=SNP",npos1) == -1) { return -1;}*/
		}
		else 
			nReturn = -1;
			
	return nReturn;
}

static int n_PickInfo2(string &str,int &nPosS,int &nPosE,Uint &urs)
{
	bool bNormal = false;
	bool bSNP = true;
	int nposEnd = 0;
	int npos1 = 0;// 自串头位置
	int npos2 = 0;// 自串尾巴位置
	int nReturn = -1;// 返回值 0:snp  1:indel -1:error
	
	
	npos2 = str.find("\t",0);
	if (npos2!= -1)
	{
		/********************************CHORM**********************************/
		  string sChr = str.substr(0,npos2);
/*		  Replace_char(sChr,"chr","");
		  if(sChr.find_first_not_of("0123456789") != -1 )
	  		return -1;
		  if(sChr == "X")
		  	nChr = 23;
		  else if(sChr == "X")
		  	nChr = 24;
		  else
		  	nChr = str2int(sChr);*/
			npos1 = npos2+1;
		/********************************POS start*************************************/
			npos2 = str.find("\t",npos1);
			string sPos= str.substr(npos1,npos2-npos1);
			int  nPos = str2int(sPos);
			nPosS = str2int(sPos);
			npos1 = npos2+1;
			nposEnd = npos1;
			/********************************POS END*************************************/
			npos2 = str.find("\t",npos1);
			string sPosE= str.substr(npos1,npos2-npos1);
			nPosE = str2int(sPosE);
			npos1= npos2+1;
		/********************************rsID*************************************/
		  string sID= "";
			npos2 = str.find("\t",npos1);
			
			sID = (npos2 == -1)?(str.substr(npos1)):(str.substr(npos1,npos2-npos1));
			
	  	if(Replace_char(sID,"rs","")!= 1)
	  		return -1;
	  	if(sID.find_first_not_of("0123456789") != -1 )
	  		return -1;
	  	urs = stoUint(sID);
	  	str = str.substr(0,nposEnd);
	  	nReturn = 0;
	  	//cout << sChr << "\t" << sPos << "\t" <<sPosE << "\t" << sID << endl;
		}
		

	return nReturn;
}

bool nGetPosInRef(Uint unRS,map<Uint,Uint>*RS_POSS,int &nPos,int &nChr)
{
	bool bFind = false;
  map<Uint,Uint>::iterator it;
	for(int i = 0;i<BIN_FILE_NUM;i++){
		nChr = i+1;
		it = RS_POSS[i].find(unRS);
		if(it != RS_POSS[i].end() )
	  {
	    nPos =  it->second;
	  	bFind = true;
	    break;
     }
  }// end for
	 
	return bFind;
	
}

bool bCheckBedErr(string sSrcfile,string sDespath,map<Uint,Uint>*RS_POSS,string &sError)
{
	bool bRe = true;
	sError = "";
	string sFilename = sSrcfile.substr(sSrcfile.rfind('/')+1,sSrcfile.length() - sSrcfile.rfind('/') -1 );
	//cout << "Check file: " << sFilename <<endl;

	string snp_file;
	snp_file = sFilename.substr(0,sFilename.rfind('.'));
	snp_file += "_error.dat";
	snp_file = sDespath + "/" +  snp_file;
	
	
	Replace_char(snp_file,"//","/");
	sError = snp_file ;
	
	cout << "Check errors file: " << snp_file <<endl;
	int  n_CountLine = 0;
	string str_Line;
	ifstream ifs;
	
	ifs.open(sSrcfile.c_str()) ;
	if (!ifs){
			cout<<"Read File Error, please check file name is right: \n"<<sSrcfile <<endl;
			sError = "读取文件出错 ";
			return false;
	}
		
	FILE *pFile;
	if ((pFile = fopen(snp_file.c_str(),"w+"))==NULL)
	{
		sError = "创建SNP.dat 文件失败 ";
		return false;
	}

	while(!ifs.eof())  
	{  
		getline(ifs,str_Line);
		//cout << str_Line << endl;
		if (0 ==str_Line.find("#")){continue;}
		//if (str_Line.find("VT=SNP") == -1) {continue;}

		n_CountLine++;
		Uint urs = 0;
		int nPosS = 0;
		int nPosE = 0;
		int nChr = -1;
		 		
 		if(0 == n_PickInfo2(str_Line,nPosS,nPosE,urs)) // 抽离出来 抽离出来位置
		{ 
			///////Get chr and pos in GRCh38
			  int nChr38 = -1;
			  int nPos38 = -1;
		    
		    if(true == nGetPosInRef(urs,RS_POSS,nPos38,nChr38) )
		    {
		      if(nPosS != nPos38 )
		    	{
		    		string sErr = str_Line + "find in reference: \tchr" + int2string(nChr38) + "\t" + int2string(nPos38) + " rs";
		    		sErr +=  Uint2string(urs) + "\n";
		    		bWritResult(pFile,&sErr); 
		    		//printf("%u \n",urs);
		    	}
		    }
		   ////////////////////////////////////////////////////////////////////////////////////////
		}					
		
	}  //end while
	
	ifs.close();  
	if (fclose (pFile) != 0)
		perror("Error occurs when close file"); 

	return bRe;
}



bool getVcfBedFile(string sSrcfile,string sDespath,string &strRe)
{
	bool bRe = true;
	strRe = "";
	string sFilename = sSrcfile.substr(sSrcfile.rfind('/')+1,sSrcfile.length() - sSrcfile.rfind('/') -1 );
	//cout << "sFilename = " << sFilename <<endl;

	string snp_file;
	snp_file = sFilename.substr(0,sFilename.rfind('.'));
	snp_file += "_SNP.bed";
	snp_file = sDespath + "/" +  snp_file;
	
	strRe = snp_file;
	Replace_char(snp_file,"//","/");
	
	//cout << "Out_put file: " << snp_file <<endl;
	int  n_CountLine = 0;
	string str_Line;
	ifstream ifs;
	
	ifs.open(sSrcfile.c_str()) ;
	if (!ifs){
			strRe = "Read File Error, please check file name is right: " +  sSrcfile;
			return false;
	}
		
	FILE *pFile;
	if ((pFile = fopen(snp_file.c_str(),"w+"))==NULL)
	{
		strRe = "Creat File Error : " +  snp_file;
		return false;
	}

	while(!ifs.eof())  
	{  
		getline(ifs,str_Line);
		if (0 ==str_Line.find("#")){continue;}
		else
		{	
			 if(0 == n_PickInfosnp(str_Line))
				{
					bWritResult(pFile,&str_Line);	
					n_CountLine++;
				}
								
		}
	}  //end while
	
	if(n_CountLine == 0)
		{
			ifs.close();// 再处理遍
			ifs.open(sSrcfile.c_str()) ;
	    if(!ifs)
	    {
			   strRe = "Read File Error, please check file name is right: " +  sSrcfile;
			   return false;
	    }
	    
	    	while(!ifs.eof())  
	        {  
		        getline(ifs,str_Line);
		        if (0 ==str_Line.find("#")){continue;}
		        else
		          {	
			         if(1 == n_PickInfosnp(str_Line)) // Normal file
				        {
					         bWritResult(pFile,&str_Line);	
					         n_CountLine++;
				        }
		         }
	        }  //end while
			
		}
		
		
		
	
	ifs.close();  
	if (fclose (pFile) != 0)
	{
			perror("Error occurs when close file"); //报告相应错误
			strRe = "Error occurs when close file:  " +  snp_file;
			return false;
	}
	
	if(n_CountLine == 0)
		{
			strRe = "[Error occurs:] Please insure there are rsID info in file: " +  sSrcfile +  "\n";
			return false;
		}
	
	


	return bRe;
}

void *  thread_fun_GRCh38BED(void* pT)  // 线程函数
{
	THREADARG  *pbg = (THREADARG  *)pT;
	//printf("thread created success ID %lu \n ",pthread_self());
	//////////////////////////////////////////////////
	string sSrcfile = pbg->sFileName;
	string str_Line;
	string strRe;
	char *cSt =(char *) strRe.data();
	
	ifstream ifs;
	ifs.open(sSrcfile.c_str());
	if (!ifs)
	{
		cout<<"Read File Error, please check file name is right: \n"<<sSrcfile <<endl;
		return (void *)cSt;
	}

	int nLinenum =0;
	Uint nPos =0;
	Uint nRs =0;
	int nChr = -1;
	while(!ifs.eof())  
	{ 
		str_Line = "";
		getline(ifs,str_Line);
		if(str_Line.length()<2)
			continue;
		if(str_Line[0] == '#')
			continue;
			
		if( 0 != n_PickRS(str_Line,nPos,nRs,nChr))
			continue;
		
		nLinenum++;

		if(nRs != 0 && nPos!=0)
			pbg->Map_PosBag->insert(make_pair(nRs,nPos));
	}//end for
	
  ifs.close();
	
	return (void *)cSt;
	/////////////////////////////////////////////////
}



int  InitRefRsPos(string strFolder,map<Uint,Uint> *mapS_Chain_RS)
{
	  ////////////////////////////////////////
	  pthread_t Init_thread[BIN_FILE_NUM];
	  THREADARG  maprsPos[BIN_FILE_NUM];
	  vector<string> aLLFiles;
    string strFileName = "";
    
	  int nCount = getAllFiles(strFolder,aLLFiles);
    if(nCount < BIN_FILE_NUM )
    	{
    		cout<< "BIN_FILE_NUM != FILE NUM";
    		return -1;
    	}
			for(int i = 0;i<nCount;i++) // creat thread
			{
				
				//ALL.chr10_GRCh38.genotypes.20170504.vcf

				strFileName = aLLFiles.at(i);
				
		    int np1 = strFileName.rfind('/')+1;
		    int np2 = strFileName.rfind('.');
		    string sfiname = strFileName.substr(np1);
		    sfiname = sfiname.substr(0,sfiname.find("_"));
		    //Replace_char(sfiname,"chr","");
		    Replace_char(sfiname,"ALL.chr","");
		    Replace_char(sfiname,"X","23");
		    Replace_char(sfiname,"Y","24");
		    if(sfiname.find_first_not_of("0123456789") != -1)
		      continue;
		    int nChr = str2int(sfiname) -1;
		    
		    //cout<< "get rs info from :" <<strFileName<< endl;
		    //cout<<nChr<< endl;
		    
		    maprsPos[nChr].nNum = nChr;
		    maprsPos[nChr].Map_PosBag = &mapS_Chain_RS[nChr];
		    maprsPos[nChr].sFileName = strFileName;
		   
			int nerror = pthread_create(&Init_thread[nChr],NULL,thread_fun_GRCh38BED,(void *)&maprsPos[nChr]);
			 if (0!=nerror)
			  	{
			  		 printf("creat thread error :    %d  %s \r\n",i,strerror(nerror));
			  		 return -1;
			  	}	
			  	
			 //sleep(10000);
			}
			
			for(int i = 0;i<BIN_FILE_NUM;i++) //wait end
			{
				char *ps = NULL;
				int nerror = pthread_join(Init_thread[i],(void **)&ps);
				if (0!=nerror)
				{
					printf("pthread_join  %s \r\n", strerror(nerror));
					return -1;
				}
				//printf("pthread_join  %s \r\n", ps);	 
			}
			//end for
	////////////////////////////////////////////////////////////////////////////////////////////
	
	cout<< "start check" << endl;
return 0;
}



void * thread_IntiMap_DAT_V2(void* pT)
	{
		THREADARGESDB *pBag = (THREADARGESDB *)pT;
		string sSrcfile = pBag->sFile;
		multimap<string,string> *mMap = pBag->pDAT_V2;
		int  n_CountLine = 0;
		string str_Line;
		string strRe = "";
		ifstream ifs;
		
		ifs.open(sSrcfile.c_str());
		if (!ifs)
		{
			//cout<<"Read File Error:\t" << sSrcfile <<endl;
		  strRe =  "Read File Error, please check file name is right: \n" + sSrcfile;
			pBag->bInit = false;
	    char *cSt =(char *) strRe.data();
	    return (void *)cSt;
		}
		while(!ifs.eof() ) 
		{  
			getline(ifs,str_Line);
			n_CountLine++;
			string skey = "";
			int npos1 = str_Line.find("\t");
			if (npos1 == -1)
				continue;
			skey = str_Line.substr(0,npos1); //1
			//cout << skey << "*"<<endl;
			mMap->insert(make_pair(skey,str_Line));
				
		}
		ifs.close();
		
	strRe =  sSrcfile + "  is ok \n";
	char *cSt =(char *) strRe.data();
	pBag->bInit = true;
	return (void *)cSt;
  }


  bool InitDATMAP_MT(string sPath) // Init format Database 
	{
		//cout<<"init \n";
		bool bRe = true;
		for (int i=0;i<CHR_NUM;i++)
		{
			string sFile = "chr" + int2string(i+1);
			sFile += ".vcf";
				
			sFile = sPath + sFile;
			
			Thrdbag[i].sFile = sFile;
			Thrdbag[i].bInit = false;
		  Thrdbag[i].pDAT_V2 = &Map_DAT_V2[i];
			
			int nerror = pthread_create(&Init_thread[i],NULL,thread_IntiMap_DAT_V2,(void *)& Thrdbag[i]);
			if (0!=nerror)
			{
				cout << "creat thread error： " << sFile << endl; 
			  printf("creat thread error :    %s \n", strerror(nerror));
			  return false;
			}
		}//end for
		
		for (int i = 0; i < CHR_NUM;i++)
			{
				char *ps = NULL;
				int nerror = pthread_join(Init_thread[i],(void **)&ps);
				if (0!=nerror)
				{
					printf("pthread_join  %s \n", strerror(nerror));
					return false;
				}
				else if( true != Thrdbag[i].bInit)
				{
					bRe = false;
					cout << "Please check file : " << Thrdbag[i].sFile << endl;
				//	printf("Iinit IMDB error  %s \n", ps);
				}
					

			}//end for
		
		return bRe;
	}
	
	
	bool bCheckFileDB(string sSrcfile,string sDespath,string &sError)
	{ 
		bool bRe = true;
	 // string sError = "";
	  string sFilename = sSrcfile.substr(sSrcfile.rfind('/')+1,sSrcfile.length() - sSrcfile.rfind('/') -1 );
	//cout << "Check file: " << sFilename <<endl;

	  string snp_file = "";
	  snp_file = sFilename.substr(0,sFilename.rfind('.'));
	  snp_file += "_db.txt";
	  snp_file = sDespath + "/" +  snp_file;
	
	  Replace_char(snp_file,"//","/");
	  cout << "Conversion errors in data base: " << snp_file <<endl;
	  int  n_CountLine = 0;
	  string str_Line;
	  ifstream ifs;
	  string strResult = "";
	
	  ifs.open(sSrcfile.c_str()) ;
	  if (!ifs){
			cout<<"Read File Error, please check file name is right: \n"<<sSrcfile <<endl;
			sError = "Read File Error: " + sSrcfile;
			return false;
	  }
		
	  FILE *pFile;
	  if((pFile = fopen(snp_file.c_str(),"w+"))==NULL){
	  	sError = "Creat File Error: " + snp_file;
		   return false;
	  }

	  while(!ifs.eof())  
	  {  
		  getline(ifs,str_Line);
		  //cout << str_Line << endl;
		  if (str_Line.length() <2){continue;}

		  n_CountLine++;
		  int npos1 = 0; 
	    int npos2 = 0; 
	    int nChr = -1;
	    string str = str_Line;
	    
	    npos1 = str.find("chr",10);
	    if (npos2 != -1)
	    {
	     npos2 = str.find("\t",npos1);
	    // cout << str << endl;
		   string sChr = str.substr(npos1,npos2 - npos1 );
		  // cout << sChr << "*" <<endl;
       Replace_char(sChr,"chr","");
		   if(sChr.find_first_not_of("0123456789XY") != -1 )
		   	{
		   		//cout << sChr << endl;
		   	  continue;
		   	}

		   if(sChr == "X")
		  	 nChr = 23;
		  else if(sChr == "Y")
		  	nChr = 24;
		  else
		  	nChr = str2int(sChr);
			 npos1 = npos2+1;
		/********************************POS start*************************************/
			npos1 = str.find(" ",npos2+1);//空格
			string sKey= str.substr(npos2+1,npos1-npos2-1);
			//cout << sKey << "*" <<endl;
			int  nPos = str2int(sKey);
			///////////////////////////////////////////////////////
			multimap<string,string>::iterator it = Map_DAT_V2[nChr-1].find(sKey);
  	  //cout << pBag->pDAT_V2[nChr-1].size()<<endl;
  	 
	     multimap<string,string>::iterator  beg,end;
	   	 if(it != Map_DAT_V2[nChr-1].end())
	   		{
	   			beg = Map_DAT_V2[nChr-1].lower_bound(sKey);
	        end = Map_DAT_V2[nChr-1].upper_bound(sKey);
	        
	        for(it = beg; it != end;it++)
	        {
	           string sRe = "";
	           if(nChr == 23)
	           	 sRe = "chrX\t";
	           else if(nChr == 24)
	           	 sRe = "chrY\t";
	           else
	           	 sRe = "chr" + int2string(nChr) + "\t";
	           	
	           sRe += it->second + "\n" ;
	           strResult += sRe;
	        }// end for 
	   		}
	   		
	   		//bWritResult(pFile,sRe);
	   		
		  }
		} // end  while
		
		//cout << strResult << endl;
		bWritResult(pFile,&strResult);
		
		 if (fclose (pFile) != 0 )
	  	{
	  		perror("Error occurs when close file"); 
	  		sError = "Error occurs when close file: ";
	  	}
		
		return true;
		 		
		
	}





int main(int arg,char *args[])
{
	 if(arg <3){
	 		cout << "Informaion: \n";
	 		cout <<sHelpInfo<<endl;
	 		return 0;
	 	}

	string sSectionDir ="";
	string strSrc = "";
	string sChainFile = "";
	string strDespath = "";
	string sRe = "";
	bool binit = false;
	bool bFolder = false;
	for(int i=0; i<arg; i++)
	{
   	string sGet = args[i];
   	switch (i)
   	{
   		case 1:  
   			     if(true ==  is_file_exist(sGet))
   			     	{
   			     	 if( true == is_dir_exist(sGet))
			 	        {
			 		        strSrc = sGet;
			 		        bFolder = true;
			 		        continue;
			 	       }
			 	       if(sGet.find(".VCF")== -1 && sGet.find(".vcf")== -1){
			 			       cout << "Source file format error < expectd format is : \".vcf\"> \n"<< endl;
			 			       cout <<sHelpInfo<<endl;
			 		         return -1;
			 		     }
				        strSrc = sGet;
			        }
			       else{ 
			       	cout << "Source file not exist " << sGet << endl; 
			       	cout <<sHelpInfo<<endl;
			       	return -1;
			       	} 
   		break;
   		case 2:
   			     if(true == is_dir_exist(sGet)) 
   			     	  sSectionDir = sGet;
			       else{ 
			       	cout << "GRCh38 vcfs folder not exist " << sGet << endl; 
			       	cout <<sHelpInfo<<endl;
			       	return -1;}
   		break;
   		case 3:
   			     if(true ==  is_file_exist(sGet) )
   			     	{
   			     		 sChainFile = sGet;
   			     		if(sChainFile.find(".chain") == -1)
  	            {
   	             cout << "stand chain file error, plz check file path\n";
   	             cout <<sHelpInfo<<endl;
   	             return 0;
  	            }
  	           
   			     	}
			       else { 
			       	cout << "chains file not exist " << sGet << endl; 
			       	cout <<sHelpInfo<<endl;
			       	return -1;}
   		break;
   		case 4:
   			{
   				 if(true == is_dir_exist(sGet)) 
   				 	 strDespath = sGet;
			      else //createDIR
			      {   
			      	if( -1 == CreateDir(sGet.c_str()))
			      	{
			      		cout << "Output folder create error " << sGet << endl; 
			      		//cout <<sHelpInfo<<endl;
			      		return -1;
			      	}
			      	strDespath = sGet;
			      	
			      }
   			}
   		break;
   		default:
   			break;	
   	}
	}
	
	if(strDespath == "")
		{
			  strDespath = "./temp/";
		    if( -1 == CreateDir(strDespath.c_str()))
			    {
			       cout << "Output folder create error " << strDespath << endl; 
			       //cout <<sHelpInfo<<endl;
			       return -1;
			    }
	      cout << "Despath :" << strDespath << endl;
		}

  

	
	/*time_t rawtime; 
  struct tm * timeinfo; 
  time ( &rawtime ); 
  timeinfo = localtime ( &rawtime ); 
  printf ( "Start time is: %d:%d:%d  \n", timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec );*/
  
  if(sChainFile.find(".chain") == -1)
  	{
   	   cout << "stand chain file error, plz check file path\n";
   	   cout <<sHelpInfo<<endl;
   	   return 0;
  	}
  	
  	map<Uint,Uint> RS_POSS[BIN_FILE_NUM]; 
  	
   string sClinvar = "./IMDB/clinvar/";
   string sGwas = "./IMDB/gwas/";
   string sHgmd = "./IMDB/hgmd/";
   string sOmim = "./IMDB/omim/";
   //string sPharmGKBDrug = "./pharmgkb_for_dat_v3/"; 
   
   if( false == InitDATMAP_MT(sClinvar) || false == InitDATMAP_MT(sGwas) || false == InitDATMAP_MT(sHgmd) ||false == InitDATMAP_MT(sOmim) ) 
   	{
   		cout<<"Initialize IMDB error, please make sure files in ./IMDB are exsit \n";
   		return -1;
   	}


  	
  	if(bFolder)
		{
			 vector<string> aLLFiles;
       string strFileName = "";
	     int nCount = getAllFiles(strSrc,aLLFiles);

			for(int i = 0;i<nCount;i++) 
			{
				strFileName = aLLFiles.at(i);
				Replace_char(strFileName,"//","/");
				if(-1 == strFileName.find(".vcf"))
					 continue;
				
				strSrc =  strFileName;
				cout << "Source file : " << strSrc << endl;
				//cout << "check file :" << strSrc << endl;
				//continue;
				if( false == getVcfBedFile(strSrc,strDespath,sRe))
   	      {
   		      cout << sRe << endl;
   		      return -1;
   	      }
        else
    	     cout << "Output bed file: " << sRe <<endl;
    	     
    	  strSrc = sRe;
        if(false == Cnvt_File_Stand(strSrc,sChainFile,strDespath,sRe))
   	     {
   		     cout << sRe << endl;
   		     return -1;
   	     }
         // 第三步 Check by GRCh38 raw vcf file
	        string strRs = "";//
	       if(binit == false)
	       	{
	       		cout << "Getting rs information from reference folder, it will takes about 30 minutes for a single chromosome and 3 hours for the human genome ...\n";
	       		if(-1 == InitRefRsPos(sSectionDir,RS_POSS))
	           {
			         cout << "Init Ref Rs Pos error\n" << endl;
			         return 0;
	           }
	          binit = true;
	       	}
	       strSrc = sRe;
	       if( false == bCheckBedErr(strSrc,strDespath,RS_POSS,sRe))
	       	{
	       		 cout << sRe << endl;
   		       return -1;
	       	}
	       	strSrc = sRe;
	       	if( false == bCheckFileDB(strSrc,strDespath,sRe))
	       		{
	       			cout << sRe << endl;
	       			continue;
	       		}
	       	cout << "\n";
		    
			}//end for
		} // folder mode
		else  // single vcf
			{
				// VCF -----> bed
				cout << "Source file : " << strSrc << endl;
        if( false == getVcfBedFile(strSrc,strDespath,sRe))
   	    {
   		    cout << sRe << endl;
   		    return -1;
   	    }
        else
    	    cout << "Output bed file: " << sRe <<endl;
        // GRCh37bed ---> GRCh38 bed
        strSrc = sRe;
        if(false == Cnvt_File_Stand(strSrc,sChainFile,strDespath,sRe))
   	     {
   		     cout << sRe << endl;
   		     return -1;
   	     }
         // Check by GRCh38 raw vcf file
         
	       cout << "Getting rs information from reference folder, it will takes about 30 minutes...\n";
	       if(-1 == InitRefRsPos(sSectionDir,RS_POSS))
	        {
			      cout << "Init error " << endl;
			      return 0;
	        }
	       strSrc = sRe;
	       //bCheckBedErr(strSrc,strDespath,RS_POSS,sRe);
	       
	       if( false == bCheckBedErr(strSrc,strDespath,RS_POSS,sRe))
	       	{
	       		 cout << sRe << endl;
   		       return -1;
	       	}
	       	strSrc = sRe;
	       	if( false == bCheckFileDB(strSrc,strDespath,sRe))
	       		{
	       			cout << sRe << endl;
	       		}  
				}
  
	return 0;
}
