function mostra_geometria(dad)
p=Scene()
for el in dad.ELEM 
  r=range(-1,stop=1,length=10)
  x = dad.NOS[el.indices,:]   # Coordenada (x,y) dos nós geométricos
  N_geo=BEM.calc_fforma.(r,Ref(el))
  ps=hcat([N_geo[i][1] for i=1:10]...)'*x
  lines!(p,ps[:,1],ps[:,2])
  if el.tipoCDC==0
    scatter!(p,x[:,1],x[:,2],color=:blue)
  elseif  el.tipoCDC==1
    scatter!(p,x[:,1],x[:,2],color=:red)
  end
end
p
end

function mostra_resultado(dad,Ts,wireframe=false)
    pts=[dad.NOS;dad.pontos_internos];
    scene=Scene()

      triin=Triangulate.TriangulateIO()
        triin.pointlist=pts'    
        (triout, vorout)=triangulate("cQ", triin)
         mesh!(scene,triout.pointlist',triout.trianglelist' , color = Ts, shading = false);
    if wireframe==true
    wireframe!(scene,scene[end][1], color = (:black, 0.6), linewidth = 3)
    end
    scene
    end