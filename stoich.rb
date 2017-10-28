#!/usr/bin/env ruby -W

require 'json'

data = ::File.open(::File.expand_path('../P_putida_model.json', __FILE__), 'r') do |f|
  source = f.read

  ::JSON.parse(source)
end

metabolite_ids = []

reaction_ids = []

stoichiometry_matrix_column_major_sparse = {}

data['metabolites'].each do |hash0|
  metabolite_id = hash0['id']

  metabolite_ids << metabolite_id
end

metabolite_ids.sort_by! { |metabolite_id|
  value = metabolite_id.split('_').last

  [value, metabolite_id]
}

data['reactions'].each do |hash0|
  reaction_id = hash0['id']
  reaction_metabolites = hash0['metabolites']
  reaction_reversibility = hash0['notes']['map_info']['reversibility']

  forward_reaction_id = '%s$f' % [reaction_id]

  reaction_ids << forward_reaction_id

  reaction_metabolites.each do |metabolite_id, stoichiometric_coefficient|
    stoichiometry_matrix_column_major_sparse[forward_reaction_id] ||= {}
    stoichiometry_matrix_column_major_sparse[forward_reaction_id][metabolite_id] = stoichiometric_coefficient
  end

  if reaction_reversibility
    backward_reaction_id = '%s$r' % [reaction_id]

    reaction_ids << backward_reaction_id

    reaction_metabolites.each do |metabolite_id, stoichiometric_coefficient|
      stoichiometry_matrix_column_major_sparse[backward_reaction_id] ||= {}
      stoichiometry_matrix_column_major_sparse[backward_reaction_id][metabolite_id] = -1 * stoichiometric_coefficient
    end
  end
end

reaction_ids.sort_by! { |reaction_id|
  value = case reaction_id.split('$').last
    when 'r' then 0
    when 'f' then 1
    else 2
  end
  
  [value, reaction_id]
}

stoichiometric_matrix_row_major_dense = ::Array.new(metabolite_ids.size) { |array| ::Array.new(reaction_ids.size, 0) }

stoichiometry_matrix_column_major_sparse.each do |reaction_id, hash0|
  reaction_column_index = reaction_ids.index(reaction_id)

  hash0.each do |metabolite_id, stoichiometric_coefficient|
    metabolite_row_index = metabolite_ids.index(metabolite_id)

    stoichiometric_matrix_row_major_dense[metabolite_row_index][reaction_column_index] = stoichiometric_coefficient
  end
end

$stdout.puts(metabolite_ids.inspect)
$stdout.puts(reaction_ids.inspect)
$stdout.puts(stoichiometric_matrix_row_major_dense.inspect)
