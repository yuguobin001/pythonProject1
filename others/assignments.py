import numpy as np
def Print_values(a,b,c):
    if a>b:
        if b>c:
            print(a,b,c)
        elif a>c:
            print(a,c,b)
        else:
            print(c,a,b)
    else :
       if b>c:
           print()
       else:
           print(c,b,a)
Print_values(np.random.randint(1,100),np.random.randint(1,100),np.random.randint(1,100))

#I got inspired by reading:https://blog.csdn.net/u011851421/article/details/83544853,https://baike.baidu.com/item/%E7%9F%A9%E9%98%B5%E4%B9%98%E6%B3%95/5446029
M1 = np.random.randint(0, 51, size=(5, 10))
M2 = np.random.randint(0, 51, size=(10, 5))
def Matrix_multip(x,y):
    m,p = x.shape
    p, n = x.shape
    result = np.zeros_like(m,n)
    for i in range(m):
        for j in range(n):
            s = 0
            for a in range(p):
                e = x[i,a]*y[a,j]
                s =s+e
            result[i,j]= s
    print(result)