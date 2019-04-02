#!/usr/bin/env gt

function usage()
  io.stderr:write("Adds a given parent attribute to all child features.\n")
  io.stderr:write(string.format("Usage: %s <attribute> < <GFF>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

function apply_to_child(root, attr, value)
  for v in root:direct_children() do
  	local a = v:get_attribute(attr)
  	if a then
  	  value = a
  	end
  	if value and not a then
  	  v:set_attribute(attr, value)
  	end
  	apply_to_child(v, attr, value)
  end
end

add_child_visitor = gt.custom_visitor_new()
function add_child_visitor:visit_feature(fn)
  apply_to_child(fn, arg[1], fn:get_attribute(arg[1]))
end

vis_stream = gt.custom_stream_new_unsorted()
function vis_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
  	node:accept(add_child_visitor)
  end
  return node
end

in_stream = gt.gff3_in_stream_new_sorted()
vis_stream.instream = in_stream
out_stream = gt.gff3_out_stream_new(vis_stream)

local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end