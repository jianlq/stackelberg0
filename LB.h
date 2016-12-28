#ifndef LB_H
#define LB_H
#include"CGraph.h"
#include <ilcplex/ilocplex.h>

inline double ILPsolve(CGraph *G, vector<demand> &eq, double &maxLinkUtil, double &sumDelay){
    IloEnv environment;
    IloModel model(environment);
    IloCplex solver(model);

    //变量
	IloArray<IloIntVarArray> x(environment, eq.size());
    for(int d = 0; d < eq.size(); d++)
		x[d] = IloIntVarArray(environment, G->m, 0, 1);
	IloNumVar z(environment, 0, IloInfinity);
	
    //优化目标
    model.add(IloMinimize(environment, z));

    //约束1
    for(int d = 0; d < eq.size(); d++)
		for(int i = 0; i < G->n; i++){
			IloExpr constraint(environment);
			for(int k = 0; k < G->adjL[i].size(); k++)
				constraint += x[d][G->adjL[i][k]->id];
			for(int k = 0; k < G->adjRL[i].size(); k++)
				constraint -= x[d][G->adjRL[i][k]->id];
			
			if(i == eq[d].org)
				model.add(constraint == 1);
			else if(i == eq[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}

    //约束2
	for(int ij = 0; ij < G->m; ij++){
		IloExpr constraint(environment);
		for(int d = 0; d < eq.size(); d++)
			constraint += eq[d].flow * x[d][ij];
		model.add(G->Link[ij]->capacity > constraint);
		model.add(z >= constraint);
	}

	solver.setOut(environment.getNullStream());
	double obj = INF;
	
    if(solver.solve()){
		obj = solver.getObjValue();
		maxLinkUtil = sumDelay = 0;
		for(int d = 0; d < eq.size(); d++)
			for(int ij = 0; ij < G->m; ij++)
				if(solver.getValue(x[d][ij]) > 0.5){
					G->Link[ij]->use += eq[d].flow;
					maxLinkUtil = max(maxLinkUtil, 1.0 * G->Link[ij]->use);
				}
			for(int ij = 0; ij < G->m; ij++)
					sumDelay += linearCal(G->Link[ij]->use,G->Link[ij]->capacity);
    }
    else{
		maxLinkUtil = sumDelay = INF;
        printf("Error: Unfeasible solve\n");
    }
	for(int i = 0; i < eq.size(); i++)
		x[i].end();
	x.end();
	environment.end();
    return obj;
}

#endif