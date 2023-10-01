# Moment-Curvature-analysis-for-UHPC-Rectangular-section


Method adopted for analysis is Fiber section mtd or Layer by layber method. In this method Entire section is divided into multiple layers and then each layer is analysed and corresponding strains, stresses, forces are calculated.
Rect beam section is as below:
                  
To calculate the Tensile stresses in the concrete, UHPC BiLinear model is considered in the analysis. Stress-Strain curve is shown ass below:
In the Bilinear analytical method, fcte and fctl are the stress points where elastic and loalization stresses are calculated using constitutive relations obtained from testing. Corresponding strains ecte and ectl are further used to determine where exactly the material is localised while applying the loads.
   
For calculation of Compressive stresses in the concrete, Unlike IS code method or conventional Hognestad’s curve, we are using UHPC stress-strain curve obtained after various testings. Graph adopted is shwn as below:
Fcm is mean targetted or max strength of UHPC. RI is the reinforcement index which is a function of fiber reinforcement.
Corresponding points are found and then used in the layer by layer analysis in the calculaton of compressive stresses in the concrete.

The basic idea is to plot the Moment v/s Curvature of the given beam section. For calculation of moment capacity we need internal forces in the UHPC and the steel reinforcement. They are function of Xu ( Neutral axis depth). 
Calculation of Xu is the only rigorous task which includes multiple iterations for attaining the force equilibrium in the section.
Cc + Cs – Ct – Ts = 0
Where,
Cc = Compressive forces in the UHPC
Cs = Compressive force because of comp steel
Ct = Tensile forces in the UHPC
Ts = Tensile force in the steel
Once we get forces, Moment is found by multiplying the forces with their corresponding lever arms.
Curvature is calculated based on the value of the Xu.
Generally heavily reinforced beams shows very less curvature indicating the less ductility.
