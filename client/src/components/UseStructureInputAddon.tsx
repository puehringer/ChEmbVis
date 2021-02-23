import * as React from "react";
import { IParticle, IParticleSelection } from "../interfaces";

export function UseStructureInputAddon({
  selection,
  setValue,
}: {
  selection?: IParticleSelection;
  setValue: (structure: string) => void;
}) {
  const selected = React.useMemo<IParticle | undefined>(
    () => Object.values(selection || {}).find((s) => s.length > 0)?.[0],
    [selection]
  );

  return (
    <div className="input-group-prepend">
      <button
        className="btn btn-outline-secondary"
        type="button"
        title="Use selected structure as input"
        disabled={!selected}
        onClick={() => (selected ? setValue(selected.structure) : undefined)}
      >
        <i className="fas fa-mouse-pointer"></i>
      </button>
    </div>
  );
}
