import ABCD
import numpy as np

class PropBeam:
	
	def __init__(self,q_in,D1,D2,D3,D4,D5,F1,F2,F3,F4):
		self.q_in = q_in
		self.D1 = D1
		self.D2 = D2
		self.D3 = D3
		self.D4 = D4
		self.D5 = D5
		self.F1 = F1
		self.F2 = F2
		self.F3 = F3
		self.F4 = F4
		self.lam = 1064e-9
		self.pi = 3.1415926

	def GetResult(self):
		q = self.q_in

		d1 = ABCD.space(self.D1)
		d2 = ABCD.space(self.D2)
		d3 = ABCD.space(self.D3)
		d4 = ABCD.space(self.D4)
		d5 = ABCD.space(self.D5)

		f1 = ABCD.lens(self.F1)
		f2 = ABCD.lens(self.F2)
		f3 = ABCD.lens(self.F3)
		f4 = ABCD.lens(self.F4) 

		total = d1 * f1 * d2 * f2 * d3 * f3 * d4 * f4 * d5
		
		q_out = (total[0,0]*q + total[0,1]) / (total[1,0]*q + total[1,1])
		psi = 1/q_out
		r_out = psi.real
		w_out = (-self.lam/(self.pi*psi.imag))**(.5)
		
		return {'q':q_out, 'r':r_out, 'w':w_out}
	
	def SpaceList(self, q_start, d):
		points = 1000		
		d1 = np.linspace(0,d,points)
		d_list = []
		q_list = []
		w_list = []
		r_list = []
		for i in range(0, points, 1):
			dz = d1[i]
			d_list.append(dz)
			abcd = ABCD.space(dz)
			q_calc = (abcd[0,0]*q_start + abcd[0,1]) / (abcd[1,0]*q_start + abcd[1,1])	
			q_list.append(q_calc)
			psi = 1/q_calc
			w_list.append((-self.lam/(self.pi*psi.imag))**(.5))
			r_list.append(1/(psi.real))
		return {'q_list':q_list, 'r_list':r_list, 'w_list':w_list, 'd_list':d_list}

	def FullModal(self):
		d_modal = []
		w_modal = []
		
		###Propgate First distance
		list1 = PropBeam.SpaceList(self,self.q_in,self.D1)
		q_list = list1['q_list']
		d_list = list1['d_list']
		w_list = list1['w_list']
		for i in range(len(d_list)):
			d_modal.append(d_list[i])
			w_modal.append(w_list[i])
		###Hits first lens
		q1 = q_list[-1]
		abcd = ABCD.lens(self.F1)
		q2 = (abcd[0,0]*q1 + abcd[0,1]) / (abcd[1,0]*q1 + abcd[1,1])
		
		####Propogate second distance
		list2 = PropBeam.SpaceList(self,q2,self.D2)
		q_list2 = list2['q_list']
		d_list2 = list2['d_list']
		w_list2 = list2['w_list']
		for i in range(len(d_list2)): 
			d_modal.append(d_list2[i]+d_list[-1])
			w_modal.append(w_list2[i])
		
		### Hits second lens
		q2 = q_list2[-1]
		abcd = ABCD.lens(self.F2)
		q3 = (abcd[0,0]*q2 + abcd[0,1]) / (abcd[1,0]*q2 + abcd[1,1])
		
		####Propogate third distance
		list3 = PropBeam.SpaceList(self,q3,self.D3)
		q_list3 = list3['q_list']
		d_list3 = list3['d_list']
		w_list3 = list3['w_list']
		for i in range(len(d_list3)): 
			d_modal.append(d_list3[i]+d_list2[-1]+d_list[-1])
			w_modal.append(w_list3[i])
		
		 ### Hits third lens
		q3 = q_list3[-1]
		abcd = ABCD.lens(self.F3)
		q4 = (abcd[0,0]*q3 + abcd[0,1]) / (abcd[1,0]*q3 + abcd[1,1])
				
		#### Propogate fourth distance
		list4 = PropBeam.SpaceList(self,q4,self.D4)
		q_list4 = list4['q_list']
		d_list4 = list4['d_list']
		w_list4 = list4['w_list']
		for i in range(len(d_list4)): 
			d_modal.append(d_list4[i]+d_list3[-1]+d_list2[-1]+d_list[-1])
			w_modal.append(w_list4[i])

		 ### Hits fourth lens
		q4 = q_list4[-1]
		abcd = ABCD.lens(self.F4)
		q5 = (abcd[0,0]*q4 + abcd[0,1]) / (abcd[1,0]*q4 + abcd[1,1])

		#### Propogate fifth distance
		list5 = PropBeam.SpaceList(self,q5,self.D5)
		q_list5 = list5['q_list']
		d_list5 = list5['d_list']
		w_list5 = list5['w_list']
		for i in range(len(d_list5)): 
			d_modal.append(d_list5[i]+d_list4[-1]+d_list3[-1]+d_list2[-1]+d_list[-1])
			w_modal.append(w_list5[i])


		return {'w_modal':w_modal, 'd_modal':d_modal}
