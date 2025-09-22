module RawInvalidations
    export raw_invalidations
    using SnoopCompileCore: @snoop_invalidations
    raw_invalidations = @snoop_invalidations using FlexUnits
end

module Shared
    module ComprehensiveToNamedTupleTrees
        export comprehensive_to_named_tuple_trees
        function comprehensive_to_named_tuple_tree_method_instance(mi)
            method = string(mi.def)
            method_instance = string(mi)
            (; method, method_instance)
        end
        function comprehensive_to_named_tuple_tree_method_instances(method_instances)
            r = []
            for mi ∈ method_instances
                t = comprehensive_to_named_tuple_tree_method_instance(mi)
                push!(r, t)
            end
            sort!(r)
        end
        function comprehensive_to_named_tuple_tree_node(node)
            method_instance = comprehensive_to_named_tuple_tree_method_instance(node.mi)
            children = comprehensive_to_named_tuple_tree_nodes(node.children)
            (; method_instance, children)
        end
        function comprehensive_to_named_tuple_tree_nodes(nodes)
            r = []
            for node ∈ nodes
                t = comprehensive_to_named_tuple_tree_node(node)
                push!(r, t)
            end
            sort!(r)
        end
        function comprehensive_to_named_tuple_tree_mt_backedges(invalidations)
            r = []
            for (type_raw, mi_or_node) ∈ invalidations
                type = string(type_raw)
                tree = if mi_or_node isa Core.MethodInstance
                    comprehensive_to_named_tuple_tree_method_instance(mi_or_node)
                else
                    comprehensive_to_named_tuple_tree_node(mi_or_node)
                end
                t = (; type, tree)
                push!(r, t)
            end
            sort!(r)
        end
        function comprehensive_to_named_tuple_tree(tree)
            method = string(tree.method)
            reason = string(tree.reason)
            mt_backedges = comprehensive_to_named_tuple_tree_mt_backedges(tree.mt_backedges)
            backedges = comprehensive_to_named_tuple_tree_nodes(tree.backedges)
            mt_cache = comprehensive_to_named_tuple_tree_method_instances(tree.mt_cache)
            mt_disable = comprehensive_to_named_tuple_tree_method_instances(tree.mt_disable)
            (; method, reason, mt_backedges, backedges, mt_cache, mt_disable)
        end
        function comprehensive_to_named_tuple_trees(trees)
            r = []
            for tree ∈ trees
                t = comprehensive_to_named_tuple_tree(tree)
                push!(r, t)
            end
            sort!(r)
        end
    end

    module ComprehensiveToNamedTuple
        export comprehensive_to_named_tuple
        using SnoopCompile: uinvalidated, invalidation_trees
        using ..ComprehensiveToNamedTupleTrees
        function comprehensive_to_named_tuple(raw_invalidations)
            invalidation_count = length(uinvalidated(raw_invalidations))
            trees = comprehensive_to_named_tuple_trees(invalidation_trees(raw_invalidations))
            (; invalidation_count, trees)
        end
    end

    module PrintJSON
        export print_json
        using JSON3: pretty
        using ..ComprehensiveToNamedTuple
        function print_json(raw_invalidations)
            pretty(comprehensive_to_named_tuple(raw_invalidations))
        end
    end

    module Exec
        using ..PrintJSON
        using ...RawInvalidations
        print_json(raw_invalidations)
    end
end