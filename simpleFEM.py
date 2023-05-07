class simpleFem:
    def __init__(self,rect):
        self.nu = 0.3 #ポアソン比
        self.E= 200*10**9 #ヤング率
        self.x=rect[0]
        self.y=rect[1]
        self.elm_num = rect[0] * rect[1] #要素数
        self.bit = np.zeros(self.elm_num) #要素のありなしのbit 0=なし 1=あり
        self.node_num = (rect[0]+1)*(rect[1]+1) #ノード数
        self.u = np.ones(self.node_num*2)  #拘束条件 1に初期化 1未知 0既知
        self.f = np.zeros(self.node_num*2) #荷重
        self.elem =  np.zeros((self.elm_num,4)) #要素の接点順番
        
        #要素の接点順番
        elem =  np.zeros((self.elm_num,4)) # 要素数×4点(四角要素)
        youso_num=0
        for i in range(self.x):#縦節点数
            for j in range(self.y): #横節点数
                self.elem[youso_num] = [i*(self.x+1)+j+1,i*(self.x+1)+j+2,i*(self.x+1)+j+self.x+3,i*(self.x+1)+j+self.x+2]
                youso_num+=1
        
    def ux(self,i,ux):
        self.u[2*i] = ux  #y方向の拘束
    
    def uy(self,i,uy):
        self.u[2*i+1] = uy  #y方向の拘束
        
    def fx(self,i,fx):
        self.f[2*i] = fx  #x方向の力
        
    def fy(self,i,fy):
        self.f[2*i+1] = fy  #y方向の力

    def calc_B(x1,y1,x2,y2,x3,y3,x4,y4,a,b):
        # D行列 平面ひずみ
        #D = np.array([[1-mu,mu,0],[mu,1-mu,0],[0,0,0.5*(1-2*mu)]])
        D = (200*10**9/(1.3*0.4)) * np.array([[0.7,0.3,0],[0.3,0.7,0],[0,0,0.2]])

        #ξ-η座標
        xi  = a/np.sqrt(3) 
        eta = b/np.sqrt(3) 

        J = (((x1-x3)*(y2-y4)-(x2-x4)*(y1-y3))+((x3-x4)*(y1-y2)-(x1-x2)*(y3-y4))*xi+((x2-x3)*(y1-y4)-(x1-x4)*(y2-y3))*eta)/8

        dN1_dx = ((+y2-y4)+(-y3+y4)*xi+(-y2+y3)*eta)/(8*J)
        dN2_dx = ((-y1+y3)+(+y3-y4)*xi+(+y1-y4)*eta)/(8*J)
        dN3_dx = ((-y2+y4)+(+y1-y2)*xi+(-y1+y4)*eta)/(8*J)
        dN4_dx = ((+y1-y3)+(-y1+y2)*xi+(+y2-y3)*eta)/(8*J)

        dN1_dy = ((-x2+x4)+(+x3-x4)*xi+(+x2-x3)*eta)/(8*J)
        dN2_dy = ((+x1-x3)+(-x3+x4)*xi+(-x1+x4)*eta)/(8*J)
        dN3_dy = ((+x2-x4)+(-x1+x2)*xi+(+x1-x4)*eta)/(8*J)
        dN4_dy = ((-x1+x3)+(+x1-x2)*xi+(-x2+x3)*eta)/(8*J)

        # B行列
        B = np.array([[dN1_dx, 0,   dN2_dx,    0,  dN3_dx,   0, dN4_dx,    0  ],
                      [   0,  dN1_dy, 0,    dN2_dy,   0,  dN3_dy,  0,   dN4_dy],
                      [dN1_dy,dN1_dx,dN2_dy,dN2_dx,dN3_dy,dN3_dx,dN4_dy,dN4_dx]])

    
        return B    

    def calc_Km(x1,y1,x2,y2,x3,y3,x4,y4):
        # D行列 平面ひずみ
        # D = np.array([[1-mu,mu,0],[mu,1-mu,0],[0,0,0.5*(1-2*mu)]])
        D = (200*10**9/(1.3*0.4)) * np.array([[0.7,0.3,0],[0.3,0.7,0],[0,0,0.2]])
        J=1
        #ξ=-1/√3 η=-1/√3
        B1 = simpleFem.calc_B(0,0,2,0,2,2,0,2,-1,-1)
        K1=B1.T@D@B1*J # K=[B.T][D][B]|J| 
    
        #ξ=1/√3 η=-1/√3
        B2 = simpleFem.calc_B(0,0,2,0,2,2,0,2,+1,-1)
        K2=B2.T@D@B2*J # K=[B.T][D][B]|J| 
    
        #ξ=-1/√3 η=1/√3
        B3 = simpleFem.calc_B(0,0,2,0,2,2,0,2,-1,+1)
        K3=B3.T@D@B3*J
    
        #ξ=1/√3 η=1/√3
        B4 = simpleFem.calc_B(0,0,2,0,2,2,0,2,1,+1)
        K4=B4.T@D@B4*J
    
        K=K1+K2+K3+K4
        return K

    # 要素マトリクス作成
    def calc_KG(self,K,ele): 
        #K = B.T @ D @ B
        KG = np.zeros((2*self.node_num,2*self.node_num))
        
        col = [ele[0],ele[0],ele[1],ele[1],ele[2],ele[2],ele[3],ele[3]] #たて
        row = [ele[0],ele[0],ele[1],ele[1],ele[2],ele[2],ele[3],ele[3]]#よこ
        idx = [0,1,0,1,0,1,0,1]
        idy = [0,1,0,1,0,1,0,1]
        a=0
        for i,y in zip(col,idy):
            b=0
            for j,x in zip(row,idx):
                KG[2*j-1+x-1,2*i-1+y-1] = K[a,b]
                b+=1
            a+=1
        return KG
    
    def calc_fem(self):
        KG = np.zeros((2*self.node_num,2*self.node_num)) #全体剛性マトリクス
        for i in range(self.elm_num):
            if fem.bit[i] == 1:
                K = np.zeros((2*4,2*4)) # 要素剛性マトリクス 
                K  += simpleFem.calc_Km(0,0,2,0,2,2,0,2) # x1,y1,x2,y2,x3,y3,x4,y4
                elm = [int(j) for j in self.elem[i]] 
                KG += simpleFem.calc_KG(self,K,elm)
            
        # 拘束条件から逆行列が存在するようにK変形
        for i in range(2*self.node_num):
            if self.u[i] == 0:
                for j in range(2*self.node_num):
                    KG[i,j] = 0
                    KG[j,i] = 0
                KG[i,i] = 1
        #逆行列
        KG_inv = np.linalg.inv(KG)

        #変位
        self.u = KG_inv @ self.f
    
    def dsp_result(self,a):    
        # 変位表示
        px=[] #変位uxを格納
        py=[] #変位uyを格納

        for i in range(self.node_num):
            px.append(self.u[2*i])
            py.append(self.u[2*i+1])
    
        #格子座標作成
        x=[] # 格子座標のx
        y=[] # 格子座標のy
        for i in range(self.y+1):
            for j in range(self.x+1):
                x.append(j)
                y.append(i)

        #グラフ表示
        #a=10**3 # 表示倍率
        for i in range(self.node_num):
            p1 = x[i]+a*px[i]
            q1 = y[i]+a*py[i]
            #plt.plot(x[i],y[i],marker='.',color = "g")
            #plt.plot(p1,q1,marker='.',color = "r")
    
        # Figureを作成
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1,1,1)
        Q = np.zeros(self.node_num)
        R = np.zeros(self.node_num)
        for i in range(self.node_num):
            Q[i]= x[i] +a*px[i]
            R[i] = y[i] +a*py[i]
        for i in range(self.elm_num):
            el = [int(i) for i in self.elem[i]]
            if self.bit[i]==1:
                b=10
                p = pat.Polygon(xy = [(Q[el[0]-1]/b, R[el[0]-1]/b), (Q[el[1]-1]/b, R[el[1]-1]/b), (Q[el[2]-1]/b, R[el[2]-1]/b),(Q[el[3]-1]/b, R[el[3]-1]/b)],fc = "lightblue", ec = "gray")
                #p = pat.Polygon(xy = [(x[el[0]-1]/50, y[el[0]-1]/50), (x[el[1]-1]/50, y[el[1]-1]/50), (x[el[2]-1]/50, y[el[2]-1]/50),(x[el[3]-1]/50, y[el[3]-1]/50)],fc = "lightblue", ec = "gray")
            
            ax.add_patch(p)
        