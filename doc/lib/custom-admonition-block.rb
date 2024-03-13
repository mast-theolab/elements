RUBY_ENGINE == 'opal' ? (require 'custom-admonition-block/extension') : (require_relative 'custom-admonition-block/extension')

Extensions.register do
  block DefinitionAdmonitionBlock
  docinfo_processor DefinitionAdmonitionBlockDocinfo
  block TheoremAdmonitionBlock
  docinfo_processor TheoremAdmonitionBlockDocinfo
  block ExerciseAdmonitionBlock
  docinfo_processor ExerciseAdmonitionBlockDocinfo
  block UsageAdmonitionBlock
  docinfo_processor UsageAdmonitionBlockDocinfo
  block ReminderAdmonitionBlock
  docinfo_processor ReminderAdmonitionBlockDocinfo
end
