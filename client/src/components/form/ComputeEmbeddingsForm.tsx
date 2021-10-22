import * as React from "react";
import { DEFAULT_COLLECTION, ICollection } from "../../interfaces";
import { embedStructures } from "../../utils/api";
import { useNameInput } from "../../utils/hooks";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { FormWrapper } from "./FormWrapper";
import { parse } from 'papaparse';

export const ComputeEmbeddingsForm = ({
  addCollection,
  loading,
  setLoading,
}: {
  addCollection(collection: ICollection): void;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [name, setName, nameInput] = useNameInput("cdddNameInput", "");
  const [cdddInput, setCdddInput] = React.useState<string>("");

  return (
    <FormWrapper
      title="Compute Embeddings"
      loading={loading}
      setLoading={setLoading}
      onSubmit={async () => {
        const structures: string[] = [];
        const additional: { [key: string]: {[key: string]: string | number | boolean} } = {};

        const result = parse<{
          smiles: string;
          [key: string]: string | number | boolean;
        }>(cdddInput, {
          header: true,
          trimHeaders: true,
          skipEmptyLines: true,
          dynamicTyping: true,
          transformHeader: (header) => {
            if(['smiles', 'structure', 'structures'].includes(header.toLocaleLowerCase())) {
              return 'smiles';
            }
            return header;
          }
        });

        if(!result.meta.fields?.includes('smiles')) {
          throw Error('No header named "smiles" found. Please include a header called "smiles" in the column of the structure.');
        }

        result.data.filter((d) => d.smiles).forEach(({smiles, ...rest}) => {
          structures.push(smiles);
          additional[smiles] = rest;
        })

        const newCollection = await embedStructures({
          structures,
          include_embedding: false
        });

        newCollection.data = newCollection.data.map((particle) => {
          const additionalProperties = additional[particle.original_structure!];
          // Inject the additional properties if there are any
          if (additionalProperties) {
            particle.properties = {
              ...(particle.properties || {}),
              ...additionalProperties,
            };
          }
          return particle;
        });
        addCollection({
          name,
          ...newCollection,
        });
        setName("");
      }}
    >
      {nameInput}
      <div className="mb-3">
        <label htmlFor="cdddTextarea">Structures to convert</label>
        <textarea
          className="form-control form-control-sm"
          id="cdddTextarea"
          rows={3}
          onChange={(e) => setCdddInput(e.currentTarget.value)}
          value={cdddInput}
        />
        <small id="cdddTextareaHelp" className="form-text text-muted">Any CSV file format with a header "smiles" is valid. You can add additional properties as columns, which will be shown in the visualizations and ranking.</small>
      </div>
      <div className="text-end">
        <ButtonWithUpload
          loading={loading}
          text="Compute"
          onUploadData={(value) => {
            setCdddInput(value || "");
            return true;
          }}
          onUploadResult={(value) => {
            if (value) {
              addCollection({
                data: JSON.parse(value),
                name: DEFAULT_COLLECTION,
              });
              return true;
            }
          }}
        />
      </div>
    </FormWrapper>
  );
};
