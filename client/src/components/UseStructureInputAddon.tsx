import * as React from "react";
import { IParticle, IParticleSelection } from "../interfaces";
import { JSMEModal } from "./JSMEModal";

export function UseStructureInputAddon({
  value,
  selection,
  setValue,
}: {
  value: string;
  selection?: IParticleSelection;
  setValue: (structure: string) => void;
}) {
  const selected = React.useMemo<IParticle | undefined>(
    () => Object.values(selection || {}).find((s) => s.length > 0)?.[0],
    [selection]
  );
  const [editorOpen, setEditorOpen] = React.useState<boolean>(false);

  return (
    <>
      <JSMEModal
        open={editorOpen}
        setOpen={setEditorOpen}
        initialSmiles={value}
        onSave={(smiles) => {
          if (smiles) {
            setValue(smiles);
          }
          setEditorOpen(false);
        }}
      />
        <button
          className="btn btn-outline-secondary"
          type="button"
          title="Use selected structure as input"
          disabled={!selected}
          onClick={() => (selected ? setValue(selected.structure) : undefined)}
        >
          <i className="fas fa-mouse-pointer"></i>
        </button>
        <button
          className="btn btn-outline-secondary"
          type="button"
          title="Draw structure"
          onClick={() => setEditorOpen(true)}
        >
          <i className="fas fa-draw-polygon"></i>
        </button>
    </>
  );
}
