require 'asciidoctor/extensions'

include Asciidoctor

# An extension that introduces a custom admonition type, complete
# with a custom icon.
#
# Usage
#
#   [QUESTION]
#   ====
#   What's the main tool for selecting colors?
#   ====
#
# or
#
#   [QUESTION]
#   What's the main tool for selecting colors?
#
class DefinitionAdmonitionBlock < Extensions::BlockProcessor
  use_dsl
  named :DEFINITION
  on_contexts :example, :paragraph

  def process parent, reader, attrs
    attrs['name'] = 'definition'
    attrs['caption'] = 'Definition'
    attrs['textlabel'] = 'Definition'
    create_block parent, :admonition, reader.lines, attrs, content_model: :compound
  end
end

class DefinitionAdmonitionBlockDocinfo < Extensions::DocinfoProcessor
  use_dsl

  def process doc
    if (doc.basebackend? 'html') && doc.backend != 'pdf'
      '<style>
.admonitionblock td.icon .icon-definition:before {content:"\f02d";color:#fff;}
</style>'
    end
  end
end

class TheoremAdmonitionBlock < Extensions::BlockProcessor
  use_dsl
  named :THEOREM
  on_contexts :example, :paragraph

  def process parent, reader, attrs
    attrs['name'] = 'theorem'
    attrs['caption'] = 'Theorem'
    attrs['textlabel'] = 'Theorem'
    create_block parent, :admonition, reader.lines, attrs, content_model: :compound
  end
end

class TheoremAdmonitionBlockDocinfo < Extensions::DocinfoProcessor
  use_dsl

  def process doc
    if (doc.basebackend? 'html') && doc.backend != 'pdf'
      '<style>
.admonitionblock td.icon .icon-theorem:before {content:"\f19d";color:#fff;}
</style>'
    end
  end
end

class ExerciseAdmonitionBlock < Extensions::BlockProcessor
  use_dsl
  named :EXERCISE
  on_contexts :example, :paragraph

  def process parent, reader, attrs
    attrs['name'] = 'exercise'
    attrs['caption'] = 'Exercise'
    attrs['textlabel'] = 'Exercise'
    create_block parent, :admonition, reader.lines, attrs, content_model: :compound
  end
end

class ExerciseAdmonitionBlockDocinfo < Extensions::DocinfoProcessor
  use_dsl

  def process doc
    if (doc.basebackend? 'html') && doc.backend != 'pdf'
      '<style>
.admonitionblock td.icon .icon-exercise:before {content:"\f059";color:#fff;}
</style>'
    end
  end
end

class UsageAdmonitionBlock < Extensions::BlockProcessor
  use_dsl
  named :USAGE
  on_contexts :example, :paragraph

  def process parent, reader, attrs
    attrs['name'] = 'usage'
    attrs['caption'] = 'Usage'
    attrs['textlabel'] = 'Usage'
    create_block parent, :admonition, reader.lines, attrs, content_model: :compound
  end
end

class UsageAdmonitionBlockDocinfo < Extensions::DocinfoProcessor
  use_dsl

  def process doc
    if (doc.basebackend? 'html') && doc.backend != 'pdf'
      '<style>
.admonitionblock td.icon .icon-usage:before {content:"\f085";color:#fff;}
</style>'
    end
  end
end

class ReminderAdmonitionBlock < Extensions::BlockProcessor
  use_dsl
  named :REMINDER
  on_contexts :example, :paragraph

  def process parent, reader, attrs
    attrs['name'] = 'reminder'
    attrs['caption'] = 'Reminder'
    attrs['textlabel'] = 'Reminder'
    create_block parent, :admonition, reader.lines, attrs, content_model: :compound
  end
end

class ReminderAdmonitionBlockDocinfo < Extensions::DocinfoProcessor
  use_dsl

  def process doc
    if (doc.basebackend? 'html') && doc.backend != 'pdf'
      '<style>
.admonitionblock td.icon .icon-reminder:before {content:"\f0f3";color:#fff;}
</style>'
    end
  end
end
