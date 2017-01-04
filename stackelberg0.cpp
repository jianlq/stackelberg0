#include"CGraph.h"
#include"LB.h"
#include"evolution.h"
#include"evolutionbit.h"

vector<demand> req;
vector<vector<demand> > pre;
bool PRE = true;
bool selfishrouting = false;


void SelfishRouting(CGraph *G, double &maxLinkUtil, double &sumDelay){
	selfishrouting = true;
	maxLinkUtil = sumDelay = 0;
	G->clearOcc();
	int NUMEVENT = req.size();

	for(int i = 0; i < NUMEVENT; i++){
		double can = G->dijkstraWeight(req[i].org, req[i].des, req[i].flow);
		if(can + SMALL >= INF){ //无解
			maxLinkUtil = INF;
			sumDelay = INF;
			return;
		}
	}
	for(int i = 0; i < G->m; i++){
		maxLinkUtil = max(maxLinkUtil, (double)G->Link[i]->use);
		sumDelay += linearCal( G->Link[i]->use,G->Link[i]->capacity);
	}
	selfishrouting = false;
}

void CentralizedOptimization(CGraph *G, double &maxLinkUtil, double &sumDelay){
	G->clearOcc();
	sumDelay = 0;
	ILPsolve(G, req, maxLinkUtil, sumDelay);
}

void LimitedRouting(CGraph *G, double &maxLinkUtil, double &sumDelay){  // 未知 + TED
	heuristicLB(G, req, req.size(), maxLinkUtil, sumDelay);

}


void StackelbergByPredict(CGraph *G, CGraph *GOR,double &maxLinkUtil, double &sumDelay, vector<demand> &predict){ //未知 + 预测 stackelberg
	PRE = true;
	int n = 100;//种群个体数目
	int m = G->m;
	evoluPopu popu(n,m,G,GOR,&req);
	evoluDiv ret = (popu).evolution();

	for(int k = 0; k < G->m; k++)
		G->Link[k]->dist = ret.x[k];

	maxLinkUtil = sumDelay = 0;
	int NUMDEMAND = req.size();
	G->clearOcc();
	for(int i = 0; i < NUMDEMAND; i++){
		double can = G->dijkstraWeight(req[i].org, req[i].des, req[i].flow);
		if(can >= INF){
			maxLinkUtil = INF;
			sumDelay = INF;
			return;
		}
	}
	for(int i = 0; i < G->m; i++){
		maxLinkUtil = max(maxLinkUtil, (double)G->Link[i]->use);
		sumDelay += linearCal( G->Link[i]->use,G->Link[i]->capacity);
	}
}

void Stackelberg(CGraph *G, CGraph *GOR,double &maxLinkUtil, double &sumDelay, vector<demand> &predict){  // 已知 + stackelberg
	PRE = false;
	int n = 100;//种群个体数目
	int m = G->m;
	evoluPopu popu(n,m,G,GOR,&req);
	evoluDiv ret = (popu).evolution();

	maxLinkUtil = ret.ability - ret.delay;
	sumDelay = ret.delay;

	if(ret.ability +SMALL >INF){
		maxLinkUtil = INF;
		sumDelay = INF;
	}

}

int TestCase(CGraph *G, CGraph *GOR,double &greedy, double &grew,double &alg1, double &alg1w,double &alg2, double &alg2w,double &alg3, double &alg3w){	

	printf("\n*****************************************\n");
	double mlu, sd;

	SelfishRouting(G, mlu, sd);  
	double dbase = sd;
	printf("SelfishRouting:\t\t%f\t%f\n", mlu, sd);
	double gremlu = mlu;

	CentralizedOptimization(G, mlu, sd);  
	double ubase = mlu;
	printf("Centralized Optimization:\t%f\t%f\n", mlu, sd);


	greedy = gremlu/ubase;
	grew = 1;

	////ilp
	alg1 = mlu/ubase;
	alg1w = sd/dbase;


	StackelbergByPredict(G, GOR,mlu, sd, req);
	printf("StackelbergByPredict:\t%f\t%f\n", mlu, sd);
	alg2 = mlu/ubase;
	alg2w = sd/dbase;

	Stackelberg(G, GOR,mlu, sd, req);
	printf("Stackelberg:\t\t%f\t%f\n", mlu, sd);

	alg3 = mlu/ubase;
	alg3w = sd/dbase;

	if(ubase >= INF)
		return 0;
	else
		return 1;
}

int main(){
	srand((unsigned)time(NULL));
	int cseNum = 200;
	int PreCase = 1;

	FILE *ret = fopen("outputFile//result.csv", "a");
	fprintf(ret, ",greedy ,,ilp,,predictS,,Stackelberg\n");
	fclose(ret);

	for(int t =0; t<5;t++){
		for(int NUMDEMAND = 5; NUMDEMAND <= 200; NUMDEMAND += 5){
			double  gremlu=0, ilpmlu = 0,spmlu = 0,spd = 0;
			double  gred = 0, ilpd = 0,smlu = 0, sd = 0;
			FILE *out = fopen("outputFile//precase.csv", "a");
			fclose(out);
			for(int numPre = 1; numPre <= PreCase ;numPre++){
				double  greedy = 0, alg1 = 0, alg2 = 0, alg3 = 0;
				double  grew = 0, alg1w = 0, alg2w = 0, alg3w = 0;
				int sucCase = 0;
				for(int cs = 0; cs < cseNum; cs++){
					int success = 0;

					genGraph(10,68,"inputFile//graph.txt");
					genGraphOR(10,5,12,"inputFile//graphOR.txt");
					CGraph *G = new CGraph("inputFile//graph.txt");
					CGraph *GOR = new CGraph("inputFile//graphOR.txt");

					//pre.clear(); //知道个数，不知道大小,猜常数
					//for(int i = 0; i < numPre; i++){
					//	req.clear();
					//	for(int i = 0; i < NUMDEMAND; i++){
					//		int s, t;
					//		do{
					//			s = rand()%G->n;
					//			t = rand()%G->n;
					//		}while(s == t || G->canNotReach(s, t));
					//		req.push_back(demand(s, t, 3));
					//	}
					//	pre.push_back(req);
					//}


					////不知道个数，不知道大小，也不知道源宿点  猜均匀分布
					pre.clear();
					req.clear();
					for(int i = 0; i < G->n; i++)
						for(int j = 0; j < G->n; j++)
							if(i != j)
								req.push_back(demand(i, j, 1));
					pre.push_back(req);


					//不知道拓扑，猜全连通图,应该是知道拓扑

					req.clear();
					for(int i = 0; i < NUMDEMAND; i++){
						int s, t;
						do{
							s = rand()%G->n;
							t = rand()%G->n;
						}while(s == t || G->canNotReach(s, t));
						req.push_back(demand(s, t, rand()%MAXDEMAND+1));
					}

					double a, b, c, d, e, f, g, h;
					success = TestCase(G,GOR, a, b, c, d, e, f, g, h);
					sucCase += success;

					if(success){
						greedy += a;
						grew += b ;
						alg1 += c;
						alg1w += d;
						alg2 += e;
						alg2w += f;
						alg3 += g;
						alg3w += h;
					}
					delete G;
					delete GOR;
				}
				out = fopen("outputFile//precase.csv", "a");
				fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f\n",numPre, greedy/sucCase,grew/sucCase,alg1/sucCase, alg1w/sucCase,alg2/sucCase, alg2w/sucCase,alg3/sucCase, alg3w/sucCase);	
				fclose(out);

				gremlu += greedy/sucCase,gred += grew/sucCase;
				ilpmlu += alg1/sucCase, ilpd += alg1w/sucCase;
				spmlu += alg2/sucCase, spd += alg2w/sucCase;
				smlu += alg3/sucCase,sd += alg3w/sucCase;
			}	

			ret = fopen("outputFile//result.csv", "a");
			fprintf(ret, "%d,%f,%f,%f,%f,%f,%f,%f,%f\n",NUMDEMAND,gremlu/PreCase, gred/PreCase, ilpmlu/PreCase, ilpd/PreCase, spmlu/PreCase,spd/PreCase,smlu/PreCase,sd/PreCase);
			fclose(ret);
		}
	}
	return 0;
}
