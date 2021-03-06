# %%

# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>

# # What is an EMU basis vector?
# 
#  [Crown and Antoniewicz, 2012](http://linkinghub.elsevier.com/retrieve/pii/S1096717611001273)  provide an algorithm for determining which fluxes are observable given a selection of tracers for 13C-metabolic flux analysis. Consider the network in figure 1.  In this system, metabolite $A$ is the substrate, metabolite $E$ is the product, and the intermediate metabolites $B, C$ and $D$ are assumed to be at metabolic and isotopic steady state. This network has two free fluxes, namely the substrate consumption rate $v_1$ and the flux of the reversible cleavage/condensation reaction $v_4$. We would like to determine flux $v_4$ based solely on the measurement of the Mass Isotope Distribution (MID) of $E$ and a known uptake rate of $A$. Three key questions that we would like to answer are: 
# 
# 1. Which labeling pattern(s) of $A$ can be used to resolve flux $v_4$? 
# 2. Are there certain tracers that will never allow flux $v_4$ to be observed? 
# 3. Is there a rational criterion that can be used a priori to identify good tracers from poor tracers? 
# 
# ## Definitions
# * **Isotopic Labeling State** represents  a substrate with a mixture of labels where the location of each labeled atom is known.   To fully specify the isotopic labeling state of the  mixture, both percentage of each component of the mixture and which atoms are labeled in each component must be described.
# * **Mass Isotopomer Distribution**  represents a substrate witha mixture of labels where the location of each labeled atom is not known. For a molecule with $n$ atoms, $M_0,M_1,\ldots,M_n$ specifies the percentage of atoms that are unlabeled, 1-labeled, up to fully labeled.
# * **Elementary Metabolite Units (EMU)** are distinct subsets of the compound's atoms.  There are 7 EMU's for a 3-atom metabolite $A$ shown in the figure from [Antoniewicz 2007](http://linkinghub.elsevier.com/retrieve/pii/S109671760600084X) below. The subscript in the first column and the shaded areas in the second column denote the atoms that are included in the EMU. The EMU size is the number of atoms included in the EMU.
# * An **EMU network** is a set of EMU's and the transitions between them.  The network can be restricted to EMUs of a particular size (Figure 3A), or it can be combined into a larger network composed of EMU's of all sizes (Figure 3B).
# *  **EMU basis vector** can be viewed as a unique way for assembling substrate EMUs to form the product. By definition, EMU basis vectors have the same EMU size as the product. Ultimately, one can interpret the Mass Isotopic Distribution (MID) of the product as a linear combination of EMU basis vector MIDs, where the coefficients are solely a function of free fluxes. The coefficients quantify the fractional contribution of each EMU basis vector to the product. Thus, by definition, the sum of the coefficients must equal one. 
# * In addition to the EMU basis vector, Crown et al also define an **EMU pathway**. Each EMU basis vector specifies a unique EMU pathway from substrate(s) to product. 
# * An **Elementary mode (ELMO)** of a network is a minimal cardinality set of reactions that form a steady-state. Minimal cardinality means that no steady-state flux distribution can be a subset of an elementary mode.  The set of elementary modes forms a basis for steady-state flux space.
# * Assertion: The elementary modes of the EMU network are the EMU pathways. Evidence for this assertion is provided below.

# %%


import cobra
dir(cobra.io.json)



import pandas as pd
import numpy as np
import scipy, sympy
from cobra import Metabolite, Reaction, Model
import cobra
import sys, re, os
sys.path.append(os.path.join(os.env['HOME'], 'Projects/src'))
from pyefm import calculate_elementary_modes, calculate_minimum_cut_sets
from d3flux import flux_map
from IPython.display import display, Image, SVG, Markdown,Latex, HTML
import sympy 
from itertools import chain

def augment_model_with_external_rxns( model, metabolites, prefix='EX_{}', sink_or_source='sink' ):
    new_model = model.copy()
    external_rxns = dict([(metabolite.id, 
                       cobra.Reaction(prefix.format(metabolite.id)))
                             for metabolite in metabolites])
    new_model.add_reactions( external_rxns.values() )
    for metabolite in metabolites:
        if sink_or_source == 'sink':
            external_rxns[metabolite.id].reaction = '{} -->'.format(metabolite.id)
        else:
            external_rxns[metabolite.id].reaction = '--> {}'.format(metabolite.id)
    return new_model
def flatten( lol ):
    return list(chain.from_iterable(lol))
def emus_of_met( met, network ):
    return flatten([m.reactions for m in network.metabolites if met in m.id])
def emus_of_product( met , network):
    return [rxn for rxn in network.reactions 
                for product in rxn.products 
                        if met in product.id]
def size_of_emu_rxn( rxn ):
    size = 0
    if len(rxn.products) == 0:
        for reactant in rxn.reactants:
            size += len(EMU(reactant.id).get_atoms()) 
    else:
        for product in rxn.products:
            size += len(EMU(product.id).get_atoms())
    return size

def extract_emu_network( old_network, size ):
    new_rxns = []
    new_network = Model('{}_size_{}'.format(old_network.id, size))
    new_network.add_reactions( [ rxn for rxn in old_network.reactions
                                        if size_of_emu_rxn( rxn ) == size ] )
    return new_network
def read_emu_network( emufilename ):
    emu_model = Model(emufilename)
    with open( emufilename ) as emu:
        emu_num = {}
        mets = set()
        for line in emu:
            rxn, eqn = line.split(',')
            rxn, eqn = rxn.strip(), eqn.strip()
            if rxn in emu_num:
                emu_num[rxn] += 1
            else:
                emu_num[rxn] = 0
            rxn_id = '{}_{}'.format( rxn, emu_num[rxn])
            emu_model.add_reactions([Reaction(rxn_id)])
            emu_model.reactions.get_by_id(rxn_id).build_reaction_from_string(eqn)
    return emu_model

class MID(object):
    """The Mass Isotopic Distribution (MID) class is just a wrapper around 
        a pandas Series with index  corresponding to 'M+0', 'M+1',...,M+{n} 
        where n is the number of atoms.  The values of the MID should sum to 100"""
    def __init__( self, mid ):
        if type(mid) is pd.core.series.Series:
            self.mid = mid
            if self.mid.sum() != 100:
                raise 
    def __getattr__( self, attr ):
        return self.mid.__getattr__( attr )
def symbolic_convolve(f, g, x, lower_limit, upper_limit):
    y = sympy.Symbol('y')
    h = g.subs(x, x - y)
    return sympyintegrate(f * h, (y, lower_limit, upper_limit))

def convolve( isotopic_labeling_state1, isotopic_labeling_state2 ):
    isotopic_labeling_state3 = np.convolve(isotopic_labeling_state1.values, isotopic_labeling_state2.values)
    idx = sorted(isotopic_labeling_state1.index.union(isotopic_labeling_state2.index))
    return pd.Series(isotopic_labeling_state3,index=idx)
    
class Tracer:
    """Tracer represents a substrate with a mixture of labels.  met is the metabolite name.  num_atoms is the total number of atoms in the metabolite.
    isotopic_labeling_state is a DataFrame with the index enumerating the number of elements of the mixture, 
    the first column is the percentage (the values of this column add up to 100), the remaining columns represent the number of carbon atoms (starting from 1). 
    The value of each element is the number of additional neutrons for the atom in each component of the mixture."""
    def __init__(self, met, num_atoms, isotopic_labeling_state ):
        self.metabolite = Metabolite(met)
        self.num_atoms = num_atoms
        self.isotopic_labeling_state = isotopic_labeling_state
        self.emu = EMU('_'.join([self.metabolite.id] + [str(i+1) for i in range(self.num_atoms )]))
    def get_metabolite( self ):
        return self.metabolite.id
    def get_isotopic_labeling_state( self ):
        return self.isotopic_labeling_state
    def get_mid( self, emu=None ):
        return self.emu.get_mid(self)
    def get_name( self, latex=True ):
        if latex:
            return ' and '.join(self.isotopic_labeling_state.index)
        
class EMU(object):
    atom_re = re.compile(r'^([A-Za-z0-9_-]+?)([0-9_]+)$')
    def __init__( self, emu ):
        self.name = emu
        self.metabolite, self.atoms, self.size = self.parse_emu(emu)
    def get_name( self,latex=True, delim=True ):
        if latex and delim:
            return '${}$'.format(self._repr_latex())
        elif latex:
            return self._repr_latex()
        else:
            return self.name
    def __str__( self ):
        return self.name
    def get_metabolite( self ):
        return self.metabolite
    def get_num_atoms( self ):
        return len( self.atoms )
    def get_atoms( self ):
        return [int(a) for a in self.atoms]
    def _repr_latex( self ):
        return r'\text{%s}_{%s}' % (self.metabolite, ','.join(self.atoms))
    def __repr__(self):
        return self._repr_latex()
    def _repr_html( self ):
        return '%s<sub>%s</sub>' % (self.metabolite, ','.join(self.atoms))
    def parse_emu( self, emu ):
        m = self.atom_re.search(emu)
        if m:
            metabolite, atoms = m.group(1), m.group(2).lstrip('_').split('_')
            size = len(atoms)
            return metabolite, atoms, size
        else:
            raise NameError('{} is not a valid EMU name that matches {}'.format(emu, self.atom_re))
    def get_isotopic_labeling_state( self, tracers ):
        return tracers[ self.get_metabolite() ].get_isotopic_labeling_state()[ self.get_atoms() ]
    def get_mid( self, tracers ):
        """tracers is a dictionary where the key is a metabolite and the value is the corresponding Tracer object.
        Each Tracer object contains the isotopic_distribution_vector, which is described above.  
        The emu contains a subset of the number of atoms in the metabolite isotopic distribution.
        By summing up the number of additional neutrons for each atom in the emu for each component of the mixture, 
        an overall Mass isotopic Distribution (MID) is obtained."""
        tracer_ist = tracers[self.get_metabolite()].get_isotopic_labeling_state()
        emu_ist = self.get_isotopic_labeling_state( tracers )
        mid = pd.Series(0, index=['M+{}'.format(i) for i in range( self.get_num_atoms() + 1 )])
        mid_idx = emu_ist[ self.get_atoms() ].sum( axis=1 ).values
        mid.iloc[ mid_idx  ] = tracers[self.get_metabolite()].get_isotopic_labeling_state()['%'].values
        mid = mid/mid.sum()*100.0
        return mid
def emus_of_rxn( rxn_name, network ):
    return ['{}: {}'.format( rxn.id, rxn.build_reaction_string()) 
          for rxn in network.reactions
          if rxn_name in rxn.id]

def sympy2array( sy ):
    return np.array(sy.tolist(), dtype=float)
def get_EMU_basis_vector_name(  emu_basis_vector ):
    """The EMU basis vector naming scheme is generated automatically from the EMU basis vector."""
    bv_name = []
    for met in emu_basis_vector.index:
        latex_met = met.split('_')
        latex_met = '%s_{%s}' % (latex_met[0], ','.join(latex_met[1:]))
        stoich = emu_basis_vector[met]
        if stoich < 0:
            if stoich == -1:
                bv_name.append( latex_met )
            else:
                bv_name.append( r'\left(' + latex_met + r'\right)^' + '{%d}' % -int(stoich))
    return '${}$'.format(r'\times '.join(bv_name))
class EMU_Basis_Vector:
    def __init__( self, emu_basis_vector ):
        self.bv = emu_basis_vector
        self.name = self.get_name(  )
    def _repr_latex( self ):
        return self.get_name().strip('$')
    def get_isotopic_labeling_state( self, tracers):
        ils = []
        for emu in self.get_input_emus(  ):
            for _ in range(np.abs(self.get_stoichiometry( emu ))):
                ils.append( emu.get_isotopic_labeling_state( tracers ))
        return ils
            
    def get_name( self,  latex=True, delim=True ):
        """The EMU basis vector naming scheme is generated automatically from the EMU basis vector."""
        bv_name = []
        
        emu_basis_vector = self.bv
        for emu in self.get_input_emus():
            emu_name = emu.get_name(latex=latex, delim=False )
            stoich = self.get_stoichiometry( emu )
            if stoich < 0:
                if stoich == -1:
                    bv_name.append( emu_name )
                else:
                    if latex:
                        bv_name.append( r'\left(' + emu_name + r'\right)^' + '{%d}' % -int(stoich))
                    else:
                        bv_name.append( r'(' + emu_name + r')**' + '%d' % -int(stoich))
        if latex and delim:
            return '${}$'.format(r'\times '.join(bv_name))
        elif latex:
            return r'\times '.join(bv_name)
        else:
            return ' x '.join( bv_name )
    def get_stoichiometry( self, emu ):
        return self.bv[emu.get_name(latex=False)].astype(int)
    def get_input_emus( self,sort_by_name_then_atoms=True ):
        if sort_by_name_then_atoms:
            return sorted([EMU( emu ) for emu in self.bv.index if self.bv[emu] < 0],reverse=True, key= lambda emu: (emu.get_metabolite(), -self.get_stoichiometry( emu )*emu.get_num_atoms()))
        else:
            return sorted([EMU( emu ) for emu in self.bv.index if self.bv[emu] < 0],reverse=True, key= lambda emu: -self.get_stoichiometry( emu )*emu.get_num_atoms())
    def get_output_emu( self ):
        return [EMU( emu ) for emu in self.bv.index if self.bv[emu > 0]]
    def convolve( self, mid1, mid2 ):
        mid3= np.convolve(mid1, mid2)/100.
        idx = ['M+{}'.format(i) for i in range(len(mid3))]
        return pd.Series(mid3,index=idx)
    def get_mid( self, tracers ):
        emus = self.get_input_emus( sort_by_name_then_atoms=False )
        # Convolve with self if stoichiometry 
        
        mid = emus[0].get_mid( tracers )
        for i in range( len( emus ) - 1 ):
            for _ in range( np.abs( self.get_stoichiometry(emus[i]) ) - 1):
                mid = self.convolve( mid, emus[i].get_mid( tracers ))
            mid = self.convolve(mid, emus[i+1].get_mid( tracers ))
        return mid    

def get_stoichiometric_matrix( network ):
    return pd.DataFrame(cobra.util.create_stoichiometric_matrix(network), 
                    index=[m.id for m in network.metabolites],
                    columns = [r.id for r in network.reactions])

def get_EMU_basis_vectors_and_pathways( network, external_rxns,external_mets, output_met, output_rxn ):
    """The EMU pathways are generated by calculating the elementary modes of the EMU Network.  The EMU basis vectors are constructed by transforming each elementary mode from a set of reactions to the corresponding set of external metabolites.  The transformation is restricted to external reactions so that the stoichiometry of the EMU basis vectors works out correctly.  The EMU pathways are normalized by the ELMO coefficient of the output reaction. The EMU basis vectors are normalized by the  EMU basis vector coefficient of the output metabolite"""
    S = get_stoichiometric_matrix( network )
    elmos = calculate_elementary_modes( network, verbose=False)
    elmos = elmos.T.divide(elmos[output_rxn]).T

    EMU_basis_vector = elmos[external_rxns].dot(S[external_rxns].T)[external_mets]
    EMU_basis_vector = EMU_basis_vector.T.divide(EMU_basis_vector[output_met]).T
    for em in EMU_basis_vector.index:
        EMU_basis_vector.loc[em,'EMU Basis Vector'] = get_EMU_basis_vector_name( EMU_basis_vector.loc[em] )
    return EMU_basis_vector, elmos
def get_nullspace( df ):
    basis_vectors = sympy.Matrix(S.values).nullspace()
    return pd.DataFrame(dict([('bv{}'.format(i), 
                        np.squeeze(np.array(bv.tolist(), dtype=float))) 
                        for i,bv in 
                            enumerate(basis_vectors)]),
                        index=df.columns)
def get_nullspace_rref( df ):
    nullspace = get_nullspace( df )
    rref, pivot_vars = sympy.Matrix(nullspace.values).rref()
    return pd.DataFrame(np.array(rref.tolist(), dtype=float),
                index= rxns), [df.columns[v] for v in pivot_vars]


# # Yeast metabolic model

emu_network = {}

emu_network['full'] = read_emu_network('yeast_VerboseFullEMUNetwork.txt')
external_mets = [met for met in emu_network['full'].metabolites 
                 if ('_e_' in met.id) and ('glc-D_e' not in met.id)]
glc_e = [met for met in emu_network['full'].metabolites 
             if 'glc-D_e' in met.id]

emu_network['full'] = augment_model_with_external_rxns(emu_network['full'],external_mets)
emu_network['full'] = augment_model_with_external_rxns(emu_network['full'],
                                                       glc_e, 
                                                       sink_or_source='source')

cobra.io.save_json_model(emu_network['full'],filename='yeast_EMU_size_{}.json'.format('full'),
                             sort=True, pretty=True)



emu_pwys = {}

size=1
emu_network[size] = extract_emu_network( emu_network['full'], size )
print("Number of EMU reactions of size {}: {}".format(size, len(emu_network[size].reactions)))
cobra.io.save_json_model(emu_network[size],filename='yeast_EMU_network_size_{}.json'.format(size),
                             sort=True, pretty=True)
emu_pwys[size] =  calculate_elementary_modes( emu_network[size],verbose=True,java_args='-Xmx80g')
emu_pwys[size].to_csv('yeast_EMU_pathways_size{}.csv'.format(size))


