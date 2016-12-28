#ifndef EVOLUTION
#define EVOLUTION
#include "Common.h"
#include"CGraph.h"

extern vector<vector<demand> > pre;
extern bool PRE;
class evoluDiv{
	private:
		static const int MUT = 5;
		static const int CUL = 30;
		static const int HOR = 20;
		CGraph *G;
		CGraph *GOR;
		vector<demand> *dem;
	public:
		vector<double> x;
		double ability;
		double delay;
		evoluDiv() {;}
		evoluDiv(int m, CGraph *g, CGraph *gor, vector<demand> *d){
			x.resize(g->m);
			G = g;
			GOR = gor;
			dem = d;
			randNature();	
		}
		evoluDiv(vector<double> &tx, CGraph *g,  CGraph *gor,vector<demand> *d){
			x.clear();
			G = g;
			GOR = gor;
			dem = d;
			for(int i = 0; i < tx.size(); i++)
				x.push_back(tx[i]);
		}

		evoluDiv mate(evoluDiv other){ 
			vector<double> nx;
			for(int i = 0; i < x.size(); i++){
				double r = 1.0 * rand() / (RAND_MAX);
				nx.push_back(r * x[i] + (1 - r) * other.x[i]);
			}
			return evoluDiv(nx, G, GOR,dem);
		}

		////ability
		void calAbility(){  
			ability = 0;
			delay = 0;
			for(int i = 0; i < G->m; i++)
				G->Link[i]->dist = x[i];

			if(PRE){
				for(int cse = 0; cse < pre.size(); cse++){
					G->clearOcc();
					double can = 0;
					for(int i = 0; i < pre[cse].size(); i++){
						double dis = G->dijkstraWeight(pre[cse][i].org, pre[cse][i].des, pre[cse][i].flow);
						can = max(can, dis);
					}

					double util = 0;
					for(int i = 0; i < G->m; i++)
						util = max(util, (double)G->Link[i]->use);
					ability += util;
					if(can + SMALL >= INF)
						ability = INF;
				}
			}
			else{
				G->clearOcc();
				double can = 0;
				for(int i = 0; i < (*dem).size(); i++){
					double dis = G->dijkstraWeight((*dem)[i].org, (*dem)[i].des, (*dem)[i].flow);
					can = max(can, dis);
				}

				double util = 0;
				for(int i = 0; i < G->m; i++){
					util = max(util, (double)G->Link[i]->use);
					delay += linearCal( G->Link[i]->use,G->Link[i]->capacity);
				}
				ability += util;
				if(can + SMALL >= INF)
					ability = INF;
				}
		}

		void randNature(){
			for(int i = 0; i < x.size(); i++){
				x[i] = 1.0 * HOR * rand() / RAND_MAX;
				if(rand()%2)
					x[i] = -x[i];
				x[i] += G->Link[i]->dist;
				x[i] = max(0.0, x[i]);
				x[i] = min(1.0*MAXWEIGHT, x[i]);
			}
		}

		void mutation(){  //变异
			for(int i = 0; i < x.size(); i++){
				x[i] += -(MAXWEIGHT * MUT/100.0) + 2 * (MAXWEIGHT * MUT/100.0) * rand() / RAND_MAX;
				x[i] = max(0.0, x[i]);
				x[i] = min(1.0 * MAXWEIGHT, x[i]);
			}
		}
		void culture(evoluDiv hero){
			for(int i = 0; i < x.size(); i++){
				double r = (CUL/100.0) * rand() / RAND_MAX;
				x[i] = (1 - r) * x[i] + r * hero.x[i];
			}
		}
};

bool evoluCmp(evoluDiv a, evoluDiv b){
	return a.ability < b.ability;
}

class evoluPopu{
	private:
		static const int YEAR = 200;
		static const int YEARCUL = 50;
		vector<evoluDiv> popu;
		CGraph *G;
		CGraph *GOR;
		vector<demand> *dem;
		double pm;
		//FILE *herofile;
	public:
		evoluDiv hero;
		// n 个体数，m：每个个体对应的解
		evoluPopu(int n, int m, CGraph *g, CGraph *gor,vector<demand> *d ){
			pm = 0.25;
			popu.clear();
			G = g;
			GOR = gor;
			dem = d;
			for(int i = 0; i < n; i++){
				evoluDiv divi(m, G, GOR,dem);
				popu.push_back(divi);
			}
			hero = evoluDiv(m, G, GOR,dem);
			//herofile = fopen("outputFile//hero.txt","w");
		}
		evoluDiv evolution(){
			//fprintf(herofile,"Start:\n ");
			/*for(int i = 0; i < hero.x.size(); i++)
				hero.x[i] = G->Link[i]->dist;*/
			
			hero.calAbility();
			//fprintf(herofile, "%f\n", hero.ability);
			//cout<< hero.ability<<endl;
			for(int i = 0; i < popu.size(); i++){
				popu[i].calAbility();
				//cout << popu[i].ability<<endl;
			}
			
			sort(popu.begin(), popu.end(), evoluCmp);
			

			for(int curYear = 1; curYear <= YEAR; curYear++){
				int n = popu.size(), getMore = 0;
				vector<evoluDiv> sons;
				for(int i = 0; i+1 < n; i+=2){
					sons.push_back(popu[i].mate(popu[i+1]));
					sons.push_back(popu[i+1].mate(popu[i]));
					sons.push_back(popu[i].mate(popu[i+1]));
				}
				int m = sons.size();
				for(int i = 0; i < m; i++){
					double p = rand()%100*0.01;
					if( p < pm ) 
						sons[i].mutation();
					if(curYear > YEARCUL)
						sons[i].culture(hero);
					sons[i].calAbility();
				}
				sort(sons.begin(), sons.end(), evoluCmp);
				popu.clear();
				for(int i = 0; i < n; i++){
					popu.push_back(sons[i]);
					if(sons[i].ability < hero.ability){
						hero = sons[i];
						getMore = 1;
					}
				}
				if(getMore){
					;
					//fprintf(herofile, "Year %d: find hero \n", curYear);
					//fprintf(herofile,"%f\n", hero.ability);
				}
			}
			//fprintf(herofile,"end\n\n\n\n");
			//fclose(herofile);
			return hero;
		}
};

#endif
