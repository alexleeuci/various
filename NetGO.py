from _collections import defaultdict
class NetGO:
    def __init__(self, g2gFile, alignFiles):
        self.pC = defaultdict(dict)
        self.CA = defaultdict(list)
        self.GOp = defaultdict(dict)
        self.pGO = defaultdict(dict)
        self.DRACONIAN = True
        self.get_pGO_GOp(g2gFile)
        for alignFile in alignFiles:
            self.get_pC_CA(alignFile)
    def SetIntersect(self, T1, T2):
        '''
        Takes pGO[protein], so a dict of GO terms, 1 meaning the protein
        has been annotated.
        Returns a dict of GO terms as well.
        '''
        out = {}
        if T1.size()>T2.size():
            for g in T2:
                if g in T1:
                    out[g]=1
        if T1.size()<=T2.size():
            for g in T1:
                if g in T2:
                    out[g]=1
        return out
    
    def K_g(self, g):
        '''
        Takes a GO term as a string, returns K(g) as int
        '''
        if(g in self.GOp):
            return len(self.GOp[g])
        else:
            return 0
    def K_gset(self, T):
        '''
        Takes pGO[p], a dict of GO terms with the names as keys, returns K(T)
        '''
        out=0
        for g in T:
            out+=self.K_g(g)
        return out
    def K_p(self, p):
        '''
        Takes a protein name as a string, returns K(p)
        '''
        if p in self.pGO:
            return self.k_gset(self.pGO[p])
        else:
            return 0
    def K_A2(self, A):
        '''
        Takes a pairwise alignment in the format A[u]=v, where u and v are pairs of proteins.
        Returns K(A)
        '''
        out = 0
        for u in A:
            if u in self.pGO:
                v = A[u]
                if v in self.pGO:
                    out+=self.K_gset(self.SetIntersect(self.pGO[u], self.pGO[v]))
        return out
    def K_AC(self, C):
        '''
        C is a dict of clusters in the format {cl:
        '''
        out = 0
        for cl in C:
            #K_C is the K(C) value for each cluster cl in C
            K_C = 0
            #M is a dict of proteins of the following format: {protein:number of times occurring in cl...}
            #T is a dict of GO terms of the following format: {GO term:number of times occurring in cl... }
            #both are across cluster cl
            M, T = defaultdict(int), defaultdict(int)
            numClusterFields=len(C[cl])
            #only calculate K_C if cluster cl has more than one protein
            if numClusterFields>1:
                u=C[cl][0]
                #if the first protein, u, is not a placeholder, then add it to M
                if u!="_" and u!="NA":
                    M[u]+=1
                    #then, for each GO term annotating protein u, add it to T
                    for g in self.pGO[u]:
                        T[g]+=1
                #Now, we iterate over all proteins in cluster cl (skipping the first one, which we already processed). Add them to M.
                for i in range(1, numClusterFields):
                    u=C[cl][i]
                    if u=="_" or u=="NA":
                        continue
                    M[u]+=1
                    #If we are in draconian mode:
                    #Check that, for each protein u in M, each GO term g in T annotates u.
                    #If the Go term g does not annotate any one of the proteins, then remove it from T.
                    if self.DRACONIAN:
                        for g in T:
                            if g not in self.pGO[u]:
                                T.pop(g)
                    #If we are not in draconian mode:
                    #Increment g's entry in T by one.
                    else:
                        for g in T:
                            self.T[g]+=1
                #If cluster cl has any annotations, and either more than one protein, or one protein that occurs more than once...
                if len(T)>0 and (len(M)>1 or (len(M)==1 and M[u]>1)):
                    if self.DRACONIAN:
                        #If we are in draconian mode, just add K_gset(T)
                        K_C+=self.K_gset(T)
                    else:
                        #If we are not in draconian mode, weed out any GO terms that annotate 1 or 0 proteins.
                        for g in T:
                            if T[g]>1:
                                K_C+=T[g]*self.K_g(g)/numClusterFields
            out+=K_C
        return out
    def sim_A2(self, A):
        return self.K_A2(A/len(self.GOp))
    def get_pC_CA(self, alignFile):
        lineCount = 0
        for cluster in alignFile:
            lineCount+=1
            for protein in cluster.split("\t"):
                #set pC
                if lineCount in self.pC[protein]:
                    self.pC[protein][lineCount]+=1
                else:
                    self.pC[protein][lineCount]=1
                #set CA
                self.CA[lineCount].append(protein)
    def get_pGO_GOp(self, g2gFile):
        for line in g2gFile:
            protein, GOterm = line.split("\t")[1:3]
            self.pGO[protein][GOterm] = 1
            if protein not in self.GOp[GOterm]:
                self.GOp[GOterm][protein]=1
            else:
                self.GOp[GOterm][protein]+=1

                






















                    